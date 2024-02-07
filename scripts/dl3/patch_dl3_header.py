from argparse import ArgumentParser
from astropy.coordinates import EarthLocation
from astropy.io import fits


 
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--dl3-file", required=True)
    args = parser.parse_args()

#    lst_location = EarthLocation.from_geodetic(-17.89139 * u.deg, 28.76139 * u.deg, 2184 * u.m)
    with fits.open(args.dl3_file, 'update') as f:
        for hdu in f:
            hdu.header['GEOLON'] = '-17.89139'
            hdu.header['GEOLAT'] = '28.76139'
            hdu.header['GEOALT'] = '2184'
