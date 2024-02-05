import logging

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import AltAz, Angle, SkyCoord
from astropy.time import Time
from gammapy.irf import Background2D, Background3D, FoVAlignment
from gammapy.maps import MapAxis, RegionGeom
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
        Apex angles of a rectangular pyramid.
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
        nbins=21,
        n_offset_bins=8,
        n_subsample=3,
        n_time_bins=5,
        offset_max="1.75 deg",
        gaussian_smoothing_3d=0.5,  # in pixels
    ):
        self.e_reco = e_reco
        self.location = location
        self.nbins = nbins
        self.n_subsample = n_subsample
        self.n_time_bins = n_time_bins
        self.n_offset_bins = n_offset_bins
        self.offset_max = Angle(offset_max)
        self.offset = MapAxis.from_bounds(
            0,
            self.offset_max,
            nbin=self.n_offset_bins,
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
        self.lon_axis_fine = MapAxis.from_bounds(
            -self.offset_max.value,
            self.offset_max.value,
            self.nbins * n_subsample,
            interp="lin",
            unit="deg",
            name="fov_lon_fine",
        )
        self.lat_axis_fine = MapAxis.from_bounds(
            -self.offset_max.value,
            self.offset_max.value,
            self.nbins * n_subsample,
            interp="lin",
            unit="deg",
            name="fov_lat_fine",
        )
        log.debug(f"Creating bkg models with e binning {self.e_reco}")
        log.debug(f"Creating bkg models with lon axis: {self.lon_axis}")
        log.debug(f"Creating bkg models with lat axis: {self.lat_axis}")
        log.debug(f"Creating bkg models with offset axis: {self.offset}")
        self.counts_map_eff = np.zeros((e_reco.nbin, nbins, nbins))
        self.counts_map_obs = np.zeros((e_reco.nbin, nbins, nbins))
        self.time_map_obs = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.time_map_eff = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.exclusion_geom = RegionGeom.from_regions(exclusion_regions)

        self.gaussian_smoothing_3d = gaussian_smoothing_3d

        log.debug("Calculating offset map from lon and lat axes")
        lon, lat = np.meshgrid(self.lon_axis.center.value, self.lat_axis.center.value)
        self.offset_map = np.sqrt(lon**2 + lat**2)

    def _fill_counts(self, obs):
        # hist events in evergy energy bin
        log.info(f"Filling counts map(s) for obs {obs.obs_id}")
        # convert coordinates from Ra/Dec to Alt/Az
        t = obs.events.time
        radec = obs.events.radec
        exclusion_mask = ~self.exclusion_geom.contains(radec)
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

    def _fill_time_maps(self, obs):
        log.info(f"Filling time map(s) for obs {obs.obs_id}")
        # time_map
        t_binning = np.linspace(obs.tstart.value, obs.tstop.value, self.n_time_bins)
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
            # Use fine axis to account for partial overlap
            az, alt = fov_to_sky(
                self.lon_axis_fine.center,
                self.lat_axis_fine.center,
                pointing_position.az,
                pointing_position.alt,
            )
            az, alt = np.meshgrid(az, alt)

            # Exclusion mask needs to be constructed for each time bin, because
            # the source will move in the FoV
            coord_radec = SkyCoord(az, alt, frame=frame).transform_to("icrs")
            # Need some reshape magic to get the correct form
            ex = ~self.exclusion_geom.contains(coord_radec)
            ex = ex.reshape(ex.shape[0], self.nbins, self.n_subsample).mean(axis=-1)
            exclusion_weight = (
                ex.T.reshape(self.nbins, self.nbins, self.n_subsample).mean(axis=-1)
            ).T

            # We have a n-deg square, so some pixels
            # are more than n-deg away from the center
            coord_lonlat = SkyCoord(lon, lat)
            mask_fov = (
                coord_lonlat.separation(SkyCoord(0 * u.deg, 0 * u.deg))
                < self.offset_max
            )
            # fill observatione time 2d arays
            observation_time_obs[mask_fov] += t_d.to(u.Unit(u.h)).value
            observation_time_eff += (
                t_d.to(u.Unit(u.h)).value * mask_fov * exclusion_weight
            )
        return u.Quantity(observation_time_eff, u.h), u.Quantity(
            observation_time_obs,
            u.h,
        )

    def _fill_all_maps(self, data_store, obs_ids=None):
        counts_obs = {}
        times_obs = {}
        counts_eff = {}
        times_eff = {}
        observations = data_store.get_observations(obs_ids, required_irf=[])
        for obs in observations:
            counts_map_eff, counts_map_obs = self._fill_counts(obs)
            time_map_eff, time_map_obs = self._fill_time_maps(obs)
            alpha_obs = time_map_eff / time_map_obs
            # remove pixels with less than half the nominal observation time
            # this avoids inflating the counts there
            threshold = 0.5
            mask_low_exposure = alpha_obs < threshold
            time_map_eff[mask_low_exposure] = 0
            counts_map_eff[:, mask_low_exposure] = 0

            counts_eff[obs.obs_id] = counts_map_eff
            counts_obs[obs.obs_id] = counts_map_obs

            times_eff[obs.obs_id] = time_map_eff
            times_obs[obs.obs_id] = time_map_obs
        return {
            "counts_obs": counts_obs,
            "counts_eff": counts_eff,
            "times_obs": times_obs,
            "times_eff": times_eff,
        }

    def run(self, data_store, obs_ids=None, cached_maps=None):
        log.info(f"running on {obs_ids}")
        observations = data_store.get_observations(obs_ids, required_irf=[])
        if cached_maps is None:
            cached_maps = self._fill_all_maps(data_store, obs_ids)
        for obs in observations:
            counts_map_eff = cached_maps["counts_eff"][obs.obs_id]
            counts_map_obs = cached_maps["counts_obs"][obs.obs_id]
            times_map_eff = cached_maps["times_eff"][obs.obs_id]
            times_map_obs = cached_maps["times_obs"][obs.obs_id]

            self.counts_map_eff += counts_map_eff
            self.counts_map_obs += counts_map_obs
            self.time_map_eff += times_map_eff
            self.time_map_obs += times_map_obs
        # 0/0 is nan...
        self.alpha_map = np.nan_to_num(
            self.time_map_eff / self.time_map_obs,
            nan=0.0,
            posinf=0,
            neginf=0,
        )
        self.bg = self._get_bg_offset()
        self.bg_rate = self.get_bg_rate()

    def _get_bg_offset_1r(self, counts_map):
        """
        For every offset bin find pixels with matching distance
        and calculate the mean bkg rate.
        """
        rmin = self.offset.edges.value[:-1] * self.offset.unit
        rmax = self.offset.edges.value[1:] * self.offset.unit
        bg_offset = []
        for rmi, rma in zip(rmin, rmax):
            mask = (self.offset_map >= rmi.value) & (self.offset_map < rma.value)
            log.debug(f"Mean over {np.count_nonzero(mask)} entries")
            sum_counts = np.sum(counts_map[mask])
            solid_angle_diff = cone_solid_angle(rma) - cone_solid_angle(rmi)

            mean_time_eff = np.mean(self.time_map_eff[mask])
            log.debug(f"Mean time: {mean_time_eff}")
            # This is actually a rate (counts/time/angle)
            counts_corrected = sum_counts / mean_time_eff / solid_angle_diff
            log.debug(f"counts: {counts_corrected} for {rmi}-{rma}")
            log.debug(f": ({solid_angle_diff}) {sum_counts}")
            bg_offset.append(counts_corrected.value)
        return np.array(bg_offset) * counts_corrected.unit

    def _get_bg_offset(self):
        """
        Create map of background rates for every bin in the time map
        """
        return [self._get_bg_offset_1r(c) for c in self.counts_map_eff]

    def get_bg_rate(self):
        """
        Divide rates by energy bin widths oto get the differential rate
        """
        bg_rate = []
        background_unit = u.Unit("s-1 MeV-1 sr-1")
        for bg_r, e_width in zip(self.bg, self.e_reco.bin_width):
            a = bg_r / e_width
            a = a.to(background_unit)
            log.debug(f"{a}")
            bg_rate.append(a)
        return bg_rate

    def get_bg_2d(self):
        # This uses self.bg_rate
        background_unit = u.Unit("s-1 MeV-1 sr-1")
        bg_2d = Background2D(
            axes=[self.e_reco, self.offset],
            data=self.bg_rate,
            unit=background_unit,
            fov_alignment=FoVAlignment.ALTAZ,
        )
        return bg_2d

    def get_bg_3d(self):
        """
        Calculate the fulld 3D bkg map.
        This works on the maps itself, not on the
        quantities calculated by the other get_bg functions
        """
        # This instead uses the counts, time and alpha maps
        bg_rate = []
        background_unit = u.Unit("s-1 MeV-1 sr-1")

        # calculate solid angle for each pixel
        # Careful: This is the bin widths this time!
        lon, lat = np.meshgrid(self.lon_axis.bin_width, self.lat_axis.bin_width)
        solid_angle_pixel = cone_solid_angle_rectangular_pyramid(lon, lat)

        bg_rate = (
            self.counts_map_eff
            / self.e_reco.bin_width.reshape(len(self.e_reco.bin_width), 1, 1)
            / solid_angle_pixel
            / self.time_map_eff
        )
        # TODO Can this be broadcasted?
        if self.gaussian_smoothing_3d:
            smoothing_kernel = Gaussian2DKernel(self.gaussian_smoothing_3d)
            for i in np.arange(len(bg_rate)):
                bg_rate[i] = convolve(bg_rate[i], smoothing_kernel)

        # nans and infs come from division by zero time, so they should be zero
        bg_rate = np.nan_to_num(bg_rate, nan=0.0, posinf=0, neginf=0)
        bg_3d = Background3D(
            axes=[self.e_reco, self.lon_axis, self.lat_axis],
            data=bg_rate.to(background_unit),
            unit=background_unit,
            fov_alignment=FoVAlignment.ALTAZ,
            meta={"FOVALIGN": "ALTAZ"},
        )
        return bg_3d
