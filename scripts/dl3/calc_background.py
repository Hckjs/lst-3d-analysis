# Code for bkgmodel adapted from Simone Mender
# https://github.com/cta-observatory/pybkgmodel/tree/add_exlcusion_region_method
import argparse
import logging
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import yaml
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord
from astropy.coordinates.erfa_astrom import ErfaAstromInterpolator, erfa_astrom
from astropy.time import Time
from gammapy.catalog import CATALOG_REGISTRY
from gammapy.data import DataStore
from gammapy.irf import Background2D, Background3D
from gammapy.maps import MapAxis, WcsGeom
from gammapy.utils.coordinates import fov_to_sky, sky_to_fov
from rich import progress

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)
erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))


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
        nbins=20,
        exclusion_radius="0.3 deg",
        offset_max="1.75 deg",
        exclusion_sources=None,
    ):
        self.e_reco = e_reco
        self.location = location
        self.nbins = nbins
        self.exclusion_radius = Angle(exclusion_radius)
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
        self.counts_map_eff = np.zeros((e_reco.nbin, nbins, nbins))
        self.counts_map_obs = np.zeros((e_reco.nbin, nbins, nbins))
        self.time_map_obs = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.time_map_eff = u.Quantity(np.zeros((nbins, nbins)), u.h)
        self.get_offset_map()
        self.exclusion_sources = exclusion_sources

    def get_offset_map(self):
        """Calculate offset to pointing position for every bin."""
        lon, lat = np.meshgrid(self.lon_axis.center.value, self.lat_axis.center.value)
        self.offset_map = np.sqrt(lon**2 + lat**2)

    def get_exclusion_mask(self, obs):
        """If no sources are defined, query 4fgl in FoV"""
        if self.exclusion_sources is not None:
            sources = self.exclusion_sources
            exclusion_mask = np.ones(len(obs.events.radec), dtype=bool)
            for s in sources:
                exclusion_mask &= s.separation(obs.events.radec) > self.exclusion_radius
        else:
            fgl = CATALOG_REGISTRY.get_cls("4fgl")()
            geom = WcsGeom.create(
                skydir=obs.pointing_radec,
                axes=[self.e_reco],
                width=2 * self.offset_max + 2 * self.exclusion_radius,
            )
            inside_geom = geom.to_image().contains(fgl.positions)
            idx = np.where(inside_geom)[0]
            exclusion_mask = (
                fgl.positions[0].separation(obs.events.radec) > self.exclusion_radius
            )
            for id in idx:
                exclusion_mask &= (
                    fgl.positions[id].separation(obs.events.radec)
                    > self.exclusion_radius
                )
        return exclusion_mask

    def fill_counts(self, obs, exclusion_mask):
        # hist events in evergy energy bin
        log.info(f"Filling counts map(s) for obs {obs.ob_id}")
        for j in range(self.e_reco.nbin):
            energy_mask = self.e_reco.edges[j] <= obs.events.energy
            energy_mask &= obs.events.energy < self.e_reco.edges[j + 1]
            mask = exclusion_mask & energy_mask
            # TODO Optimize mask usage
            # You could avoid some transformations here
            #
            # convert coordinates from Ra/Dec to Alt/Az
            t = obs.events.time
            frame = AltAz(obstime=t, location=self.location)
            pointing_altaz = obs.events.pointing_radec.transform_to(frame)
            position_events = obs.events.radec.transform_to(frame)
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
            #
            if j == 0:
                counts_map_eff = counts_eff
                counts_map_obs = counts_obs
            else:
                counts_map_eff = np.dstack((counts_map_eff, counts_eff))
                counts_map_obs = np.dstack((counts_map_obs, counts_obs))
        counts_map_eff = counts_map_eff.transpose()
        counts_map_obs = counts_map_obs.transpose()
        self.counts_map_eff += counts_map_eff
        self.counts_map_obs += counts_map_obs

    def fill_time_maps(self, obs):
        log.info(f"Filling time map(s) for obs {obs.ob_id}")
        # define exclusion mask for all sources in 4fgl catalog in the region
        fgl = CATALOG_REGISTRY.get_cls("4fgl")()
        geom = WcsGeom.create(
            skydir=obs.pointing_radec,
            axes=[self.e_reco],
            width=2 * self.offset_max + 2 * self.exclusion_radius,
        )
        inside_geom = geom.to_image().contains(fgl.positions)
        idx = np.where(inside_geom)[0]

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

            coord_radec = SkyCoord(az, alt, frame=frame).transform_to("icrs")
            # calculate masks for FoV and exclusion regions in Ra/Dec coordinates
            coord_lonlat = SkyCoord(lon, lat)
            mask_fov = (
                coord_lonlat.separation(SkyCoord(0 * u.deg, 0 * u.deg))
                < self.offset_max
            )
            exclusion_mask = (
                fgl.positions[0].separation(coord_radec) > self.exclusion_radius
            )
            for id in idx:
                exclusion_mask &= (
                    fgl.positions[id].separation(coord_radec) > self.exclusion_radius
                )
            # fill observatione time 2d arays
            observation_time_obs[mask_fov] += t_d.to(u.Unit(u.h)).value
            observation_time_eff[exclusion_mask & mask_fov] += t_d.to(u.Unit(u.h)).value
        self.time_map_obs += u.Quantity(observation_time_obs, u.h)
        self.time_map_eff += u.Quantity(observation_time_eff, u.h)

    def run(self, data_store, obs_ids=None):
        observations = data_store.get_observations(obs_ids, required_irf=[])
        log.info("Creating background map")
        for obs in progress.track(observations):
            exclusion_mask = self.get_exclusion_mask(obs)
            self.fill_counts(obs, exclusion_mask)
            self.fill_time_maps(obs)
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


def main():
    """
    Function running the entire background reconstruction procedure.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", default="bkg")
    parser.add_argument("--dummy-output", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    out = Path(args.output_dir)

    with open(args.config) as f:
        config = yaml.safe_load(f)
    e_binning = config["binning"]["energy"]
    fov_binning = config["binning"]["offset"]
    exclusion = config["exclusion"]

    # TODO Define that properly somewhere
    location = EarthLocation.of_site("Roque de los Muchachos")
    e_reco = MapAxis.from_energy_bounds(
        u.Quantity(e_binning["min"]),
        u.Quantity(e_binning["ax"]),
        e_binning["n_bins"],
        name="energy",
    )
    bkg_maker = ExclusionMapBackgroundMaker(
        e_reco,
        location,
        nbins=fov_binning["n_bins"],
        offset_max=u.Quantity(fov_binning["max"]),
        exclusion_radius=u.Quantity(exclusion["radius"]),
        exclusion_sources=[SkyCoord(**s) for s in exclusion["sources"]],
    )

    ds = DataStore.from_dir(args.input_dir)

    # Select similar runs. This is only zenith right now
    # Ra/dec is assumed to be correclty incorporated earlier
    # Az is neglected for now, maybe thats fine for one LST
    # Time is neglected as well. Very different runs should be excluded before this
    # and there is no notion of MC-periods in LST so far
    criteria = pd.DataFrame(
        {
            "obs_id": ds.obs_table["OBS_ID"],
            "cos_zenith": np.cos(ds.obs_table["ZEN_PNT"].quantity),
        },
    )
    for obs_id in ds.obs_ids:
        cos_zenith_diff = np.abs(
            (criteria["cos_zenith"] - criteria["cos_zenith"][0]).values,
        )
        mask = cos_zenith_diff < config["max_cos_zenith_diff"]
        selected_ids = criteria["obs_id"][mask].data
        # select fitting runs
        bkg_maker.run(ds, selected_ids)
        if config["hdu_type"] == "3D":
            bkg = bkg_maker.get_bg_3d()
        elif config["hdu_type"] == "2D":
            bkg = bkg_maker.get_bg_2d()
        else:
            raise NotImplementedError()
        bkg.write(
            out / f"{config['prefix']}_{obs_id}.fits.gz",
            overwrite=args.overwrite,
        )
    Path(args.dummy_output).touch()


if __name__ == "__main__":
    main()
