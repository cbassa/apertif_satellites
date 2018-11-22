#!/usr/bin/env python
from __future__ import print_function
import argparse
import datetime
import ephem
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import numpy as np

# python 2 vs 3 importing
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen


class two_line_element:
    """TLE class"""

    def __init__(self, tle0, tle1, tle2):
        """Define a tle"""

        self.tle0 = tle0.decode("utf-8").strip()
        self.tle1 = tle1.decode("utf-8")
        self.tle2 = tle2.decode("utf-8")
        if self.tle0[:2] == "0 ":
            self.name = self.tle0[2:].strip()
        else:
            self.name = self.tle0.strip()
        self.id = self.tle1.split(" ")[1][:5]


def download_tles(urlstr="https://celestrak.com/NORAD/elements/galileo.txt"):
    """Download two-line elements from the specified url."""
    """Returns a list of two_line_element objects."""
    lines = urlopen(urlstr).readlines()
    return [two_line_element(lines[i], lines[i+1], lines[i+2])
            for i in range(0, len(lines), 3)]


if __name__ == "__main__":
    # Read command line arguments
    parser = argparse.ArgumentParser(
        description="Compute Galileo satellite positions for WSRT")
    parser.add_argument("-t", "--time",
                        help="Time (YYYY-MM-DD HH:MM:SS) [default: now]",
                        default=datetime.datetime.utcnow().strftime(
                            "%Y-%m-%dT%H:%M:%S"
                        ))
    args = parser.parse_args()

    # Set start and end times
    tobs = datetime.datetime.strptime(args.time, "%Y-%m-%dT%H:%M:%S")

    # Download TLEs
    tles = download_tles()

    # Set observer
    observer = ephem.Observer()
    observer.lon = "6.60334"
    observer.lat = "52.91474"
    observer.elevation = 0
    location = EarthLocation(lat=52.91474*u.deg, lon=6.60334*u.deg, height=0*u.m)

    print("Satellites are The Worst")
    print("Positions of Galileo satellites as seen by WSRT at", tobs)
    print("\n==================================================\n")
    print("Satellite          | R.A.         |  Decl.       |  Alt. (deg)")
    
    # Loop over TLEs
    for tle in tles:
        # Set satellite
        satellite = ephem.readtle(str(tle.tle0),
                                  str(tle.tle1),
                                  str(tle.tle2))
        
        # Set observer date
        observer.date = ephem.date(tobs)

        # Compute satellite observer
        satellite.compute(observer)

        # Convert to skycoord
        pos = SkyCoord(ra=float(satellite.ra), dec=float(satellite.dec), unit="rad", frame="icrs")
        pos_altaz = pos.transform_to(AltAz(obstime=tobs, location=location))

        if pos_altaz.alt.degree>0.0:
            print("%s %30s %8.1f %s"%(tle.tle0, pos.to_string("hmsdms"), pos_altaz.alt.degree, tle.id))
