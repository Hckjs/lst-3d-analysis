import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, angular_separation
from astropy.time import Time


def sin_delta(altaz: AltAz):
    """Delta is the angle between pointing and magnetic field."""
    # Values from
    # https://geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html
    # for La Palma coordinates and date=2021-12-01.
    #
    # Also used here:
    # https://github.com/cta-observatory/lst-sim-config/issues/2
    b_inc = u.Quantity(-37.36, u.deg)
    b_dec = u.Quantity(-4.84, u.deg)

    delta = angular_separation(
        lon1=altaz.az,
        lat1=altaz.alt,
        lon2=b_dec,
        lat2=b_inc,
    )
    return np.sin(delta.to_value(u.rad))


def cos_zenith(altaz: AltAz):
    return np.cos(altaz.zen.to_value(u.rad))


@u.quantity_input
def build_altaz(*, alt: u.deg = None, zd: u.deg = None, az: u.deg = None) -> AltAz:
    """Build AltAz from zenith distance and azimuth.

    This function exists to make no mistakes when translating
    zenith distance to altitude.

    Altitude takes precedence over zenith, if given both.

    location and obstime is needed for transformations with
    astropy version 4.3.1

    Obstime is fixed.
    """
    if alt is None:
        if zd is None:
            raise ValueError("Specify either alt or zd")
        zenith = u.Quantity(90, u.deg)
        alt = zenith - zd

    location = EarthLocation.from_geodetic(
        u.Quantity(-17.89139, u.deg),
        u.Quantity(28.76139, u.deg),
        height=u.Quantity(2184, u.m),
    )
    obstime = Time("2022-01-01T00:00")

    return AltAz(alt=alt, az=az, location=location, obstime=obstime)


def get_theta_az_from_node(node: str) -> np.ndarray:
    """Strip theta and az from name of node.

    `node` must be of kind 'node_theta_10.0_az_102.199_'.
    """
    _, _, theta, _, az, _ = node.split("_")
    return np.array([theta, az], dtype=float)


def get_pointings_of_irfs(filelist) -> AltAz:
    """From the list of directory names with AllSky IRFs,
    build the AltAz frame with pointings.

    The names are of kind 'node_theta_10.0_az_102.199_'.
    """
    theta, az = np.array([get_theta_az_from_node(f) for f in filelist]).T

    return build_altaz(zd=theta * u.deg, az=az * u.deg)


def euclidean_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
