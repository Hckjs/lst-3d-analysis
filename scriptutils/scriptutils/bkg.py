import logging

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, Angle, SkyCoord
from astropy.time import Time
from gammapy.irf import Background2D, Background3D
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.utils.coordinates import fov_to_sky, sky_to_fov

log = logging.getLogger(__name__)


def cone_solid_angle(angle):
    """
    Calculate the solid angle of a view cone.
    Parameters
    ----------
    angle: astropy.units.Quantity or astropy.coordinates.Angle
        Opening angle of the view cone.
    Returns
    -------
    solid_angle: astropy.units.Quantity
        Solid angle of a view cone with opening angle ``angle``.
    """
    angle = angle.to(u.Unit("rad"))
    solid_angle = 2 * np.pi * (1 - np.cos(angle)) * u.sr
    return solid_angle


def cone_solid_angle_rectangular_pyramid(a, b):
    """
    Calculate the solid angle of a view cone.
    Parameters
    ----------
    angle: astropy.units.Quantity or astropy.coordinates.Angle
        Opening angle of the view cone.
    Returns
    -------
    solid_angle: astropy.units.Quantity
        Apex angles of a retangluar pyramid.
    """
    a = a.to(u.Unit("rad"))
    b = b.to(u.Unit("rad"))
    solid_angle = 4 * np.arcsin(np.sin(a / 2) * np.sin(b / 2)).value * u.sr
    return solid_angle


class ExclusionMapBackgroundMaker:
    """Exclusion map background algorithm.
    Calculates background in FOV coordinate system aligned with the `ALTAZ` system.
    """

    def __init__(  # noqa: PLR0913
        self,
        e_reco,
        location,
        exclusion_regions,
        nbins=20,
        offset_max="1.75 deg",
    ):
        self.e_reco = e_reco
        self.location = location
        self.nbins = nbins
        self.offset_max = Angle(offset_max)
        self.offset = MapAxis.from_bounds(
            0,
            self.offset_max,
            nbin=nbins,
            interp="lin",
            unit="deg",
            name="offset",
        )
        self.lon_axis = MapAxis.from_bounds(
            -self.offset_max.value,
            self.offset_max.value,
            self.nbins,
            interp="lin",
            unit="deg",
            name="fov_lon",
        )
        self.lat_axis = MapAxis.from_bounds(
            -self.offset_max.value,
            self.offset_max.value,
            self.nbins,
            interp="lin",
            unit="deg",
            name="fov_lat",
        )

        log.debug(f"Creating bkg models with e binning {self.e_reco}")
        log.debug(f"Creating bkg models with lon axis: {self.lon_axis}")
        log.debug(f"Creating bkg models with lat axis: {self.lat_axis}")
        self.counts_map_eff = np.zeros((e_reco.nbin, nbins, nbins))
        self.counts_map_obs = np.zeros((e_reco.nbin, nbins, nbins))
        self.time_map_obs = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.time_map_eff = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.get_offset_map()
        self.exclusion_geom = RegionGeom.from_regions(exclusion_regions)

    def get_offset_map(self):
        """Calculate offset to pointing position for every bin."""
        lon, lat = np.meshgrid(self.lon_axis.center.value, self.lat_axis.center.value)
        self.offset_map = np.sqrt(lon**2 + lat**2)

    def get_exclusion_mask(self, obs):
        radec = obs.events.radec
        # select all
        mask = np.ones(len(radec), dtype=bool)
        if self.exclusion_mask:
            mask &= ~self.exclusion_geom.contains(radec)
        return mask

    def fill_counts(self, obs, exclusion_mask):
        # hist events in evergy energy bin
        log.info(f"Filling counts map(s) for obs {obs.obs_id}")
        # convert coordinates from Ra/Dec to Alt/Az
        t = obs.events.time
        frame = AltAz(obstime=t, location=self.location)
        pointing_altaz = obs.events.pointing_radec.transform_to(frame)
        position_events = obs.events.radec.transform_to(frame)
        for j in range(self.e_reco.nbin):
            energy_mask = self.e_reco.edges[j] <= obs.events.energy
            energy_mask &= obs.events.energy < self.e_reco.edges[j + 1]
            mask = exclusion_mask & energy_mask
            log.debug(
                f"Transforming {np.count_nonzero(energy_mask)} events in bin {j}",
            )
            # convert Alt/Az to Alt/Az FoV
            # effective counts
            lon, lat = sky_to_fov(
                position_events.az[mask],
                position_events.alt[mask],
                pointing_altaz.az[mask],
                pointing_altaz.alt[mask],
            )
            counts_eff, xedges, yedges = np.histogram2d(
                lon.value,
                lat.value,
                bins=(self.lon_axis.edges.value, self.lat_axis.edges.value),
            )
            # observed counts
            lon, lat = sky_to_fov(
                position_events.az[energy_mask],
                position_events.alt[energy_mask],
                pointing_altaz.az[energy_mask],
                pointing_altaz.alt[energy_mask],
            )
            counts_obs, xedges, yedges = np.histogram2d(
                lon.value,
                lat.value,
                bins=(self.lon_axis.edges.value, self.lat_axis.edges.value),
            )
            if j == 0:
                counts_map_eff = counts_eff
                counts_map_obs = counts_obs
            else:
                counts_map_eff = np.dstack((counts_map_eff, counts_eff))
                counts_map_obs = np.dstack((counts_map_obs, counts_obs))
        counts_map_eff = counts_map_eff.transpose()
        counts_map_obs = counts_map_obs.transpose()
        return counts_map_eff, counts_map_obs

    def fill_time_maps(self, obs, exclusion_mask):
        log.info(f"Filling time map(s) for obs {obs.obs_id}")
        WcsGeom.create(
            skydir=obs.pointing_radec,
            axes=[self.e_reco],
            width=2 * self.offset_max + 2 * self.exclusion_radius,
        )

        # time_map
        t_binning = np.linspace(obs.tstart.value, obs.tstop.value, 30)
        t_binning = Time(t_binning, format="mjd")
        t_delta = np.diff(t_binning)
        t_center = t_binning[:-1] + 0.5 * t_delta
        lon, lat = np.meshgrid(self.lon_axis.center, self.lat_axis.center)
        # create observation time 2d arrays
        observation_time_obs = np.zeros((self.nbins, self.nbins))
        observation_time_eff = np.zeros((self.nbins, self.nbins))
        # iterate time bins
        pointing_positions = []
        for t_c, t_d in zip(t_center, t_delta):
            # transform from camera coordinates to FoV coordinates (Alt/Az)
            # dependent on the time
            frame = AltAz(obstime=t_c, location=self.location)
            pointing_position = obs.pointing_radec.transform_to(frame)
            pointing_positions.append(pointing_position)
            az, alt = fov_to_sky(
                self.lon_axis.center,
                self.lat_axis.center,
                pointing_position.az,
                pointing_position.alt,
            )
            az, alt = np.meshgrid(az, alt)

            SkyCoord(az, alt, frame=frame).transform_to("icrs")
            # calculate masks for FoV and exclusion regions in Ra/Dec coordinates
            coord_lonlat = SkyCoord(lon, lat)
            mask_fov = (
                coord_lonlat.separation(SkyCoord(0 * u.deg, 0 * u.deg))
                < self.offset_max
            )
            # fill observatione time 2d arays
            observation_time_obs[mask_fov] += t_d.to(u.Unit(u.h)).value
            observation_time_eff[exclusion_mask & mask_fov] += t_d.to(u.Unit(u.h)).value
        return u.Quantity(observation_time_eff, u.h), u.Quantity(
            observation_time_obs,
            u.h,
        )

    def fill_all_maps(self, data_store, obs_ids=None):
        counts_obs = {}
        times_obs = {}
        counts_eff = {}
        times_eff = {}
        observations = data_store.get_observations(obs_ids, required_irf=[])
        for obs in observations:
            exclusion_mask = self.get_exclusion_mask(obs)

            counts_map_eff, counts_map_obs = self.fill_counts(obs, exclusion_mask)
            counts_eff[obs.obs_id] = counts_map_eff
            counts_obs[obs.obs_id] = counts_map_obs

            time_map_eff, time_map_obs = self.fill_time_maps(obs, exclusion_mask)
            times_eff[obs.obs_id] = time_map_eff
            times_obs[obs.obs_id] = time_map_obs
        return {
            "counts_obs": counts_obs,
            "counts_eff": counts_eff,
            "times_obs": times_obs,
            "times_eff": times_eff,
        }

    def run(self, data_store, obs_ids=None, cached_maps=None):
        print(f"running on {obs_ids}")
        observations = data_store.get_observations(obs_ids, required_irf=[])
        if cached_maps is None:
            cached_maps = self.fill_all_maps(data_store, obs_ids)
        for obs in observations:
            counts_map_eff = cached_maps["counts_eff"][obs.obs_id]
            counts_map_obs = cached_maps["counts_obs"][obs.obs_id]
            times_map_eff = cached_maps["times_eff"][obs.obs_id]
            times_map_obs = cached_maps["times_obs"][obs.obs_id]

            self.counts_map_eff += counts_map_eff
            self.counts_map_obs += counts_map_obs
            self.time_map_eff += times_map_eff
            self.time_map_obs += times_map_obs
        self.alpha_map = self.time_map_eff / self.time_map_obs
        self.bg = self.get_bg_offset(self.counts_map_eff)
        self.bg_rate = self.get_bg_rate()

    def get_bg_offset_1r(self, counts_map):
        rmin = self.offset.edges.value[:-1] * self.offset.unit
        rmax = self.offset.edges.value[1:] * self.offset.unit
        bg_offset = []
        for rmi, rma in zip(rmin, rmax):
            mask = (self.offset_map >= rmi.value) & (self.offset_map < rma.value)
            sum_counts = np.sum(counts_map[mask])
            solid_angle_diff = cone_solid_angle(rma) - cone_solid_angle(rmi)
            mean_alpha = np.mean(self.alpha_map[mask])
            mean_time = np.mean(self.time_map_obs)
            counts_corrected = sum_counts / mean_alpha / solid_angle_diff / mean_time
            bg_offset.append(counts_corrected.value)
        return np.array(bg_offset) * counts_corrected.unit

    def get_bg_offset(self, counts_map):
        return [self.get_bg_offset_1r(c) for c in self.counts_map_eff]

    def get_bg_rate(self):
        bg_rate = []
        background_unit = u.Unit("s-1 MeV-1 sr-1")
        for bg_r, e_width in zip(self.bg, self.e_reco.bin_width):
            a = bg_r / e_width
            a = a.to(background_unit)
            bg_rate.append(a)
        return bg_rate

    def get_bg_2d(self):
        background_unit = u.Unit("s-1 MeV-1 sr-1")
        bg_2d = Background2D(
            axes=[self.e_reco, self.offset],
            data=self.bg_rate,
            unit=background_unit,
        )
        return bg_2d

    def get_bg_3d(self):
        bg_3d_counts = self.counts_map_eff / self.alpha_map
        bg_rate = []
        background_unit = u.Unit("s-1 MeV-1 sr-1")
        # calculate solid angle for each pixel
        lon, lat = np.meshgrid(self.lon_axis.bin_width, self.lat_axis.bin_width)
        solid_angle_pixel = cone_solid_angle_rectangular_pyramid(lon, lat)
        # go through every energy bin
        for bg_r, e_width in zip(bg_3d_counts, self.e_reco.bin_width):
            a = bg_r / e_width / self.time_map_obs / solid_angle_pixel
            a[np.isnan(a)] = 0
            a = a.to(background_unit)
            bg_rate.append(a)
        bg_3d = Background3D(
            axes=[self.e_reco, self.lon_axis, self.lat_axis],
            data=bg_rate,
            unit=background_unit,
            meta={"FOVALIGN": "ALTAZ"},
        )
        return bg_3d
