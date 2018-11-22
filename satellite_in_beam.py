#!/usr/bin/env python
from __future__ import print_function
import argparse
import datetime
import ephem
from astropy.coordinates import SkyCoord
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

        self.tle0 = tle0.decode("utf-8")
        self.tle1 = tle1.decode("utf-8")
        self.tle2 = tle2.decode("utf-8")
        if self.tle0[:2]=="0 ":
            self.name = self.tle0[2:].strip()
        else:
            self.name = self.tle0.strip()
        self.id = self.tle1.split(" ")[1][:5]


def download_tles(urlstr="https://celestrak.com/NORAD/elements/beidou.txt"):
    """Download two-line elements from the specified url."""
    """Returns a list of two_line_element objects."""
    lines = urlopen(urlstr).readlines()
    return [two_line_element(lines[i], lines[i+1], lines[i+2]) for i in range(0, len(lines), 3)]

if __name__ == "__main__":
    # Read command line arguments
    parser = argparse.ArgumentParser(description="Compute if a satellite will pass through the telescope beam")
    parser.add_argument("-t", "--starttime",
                        help="Start time (YYYY-MM-DD HH:MM:SS) [default: now]",
                        default=datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"))
    parser.add_argument("-l", "--length", help="Search length (seconds) [default: 600s]", type=int, default=600)
    parser.add_argument("-R", "--ra", help="Right ascension (deg)", type=float)
    parser.add_argument("-D", "--dec", help="Declination (deg)", type=float)
    parser.add_argument("-r", "--radius", help="Search radius (deg) [default: 1 deg]", type=float, default=1)
    parser.add_argument("-d", "--dt", help="Time step (seconds) [default: 10s]", type=int, default=10)
    args = parser.parse_args()

    # Set start and end times
    tstart = datetime.datetime.strptime(args.starttime, "%Y-%m-%dT%H:%M:%S")
    tend = tstart+datetime.timedelta(seconds=args.length)

    # Time range
    dt = range(0, args.length, args.dt)
    
    # Set position
    beam = SkyCoord(ra=args.ra, dec=args.dec, unit="deg", frame="icrs")

    # Download TLEs
    tles = download_tles()
    
    # Set observer
    observer = ephem.Observer()
    observer.lon = "6.60334"
    observer.lat = "52.91474"
    observer.elevation = 0

    # Loop over TLEs
    for tle in tles:
        # Set satellite
        satellite = ephem.readtle(str(tle.tle0),
                                  str(tle.tle1),
                                  str(tle.tle2))
        
        # Loop over time
        ra = []
        dec = []
        for t in dt:
            # Set observer date
            observer.date = ephem.date(tstart+datetime.timedelta(seconds=t))

            # Compute satellite observer
            satellite.compute(observer)
            ra.append(float(satellite.ra))
            dec.append(float(satellite.dec))

        # Convert to skycoord
        pos = SkyCoord(ra=ra, dec=dec, unit="rad", frame="icrs")

        # Offsets
        r = beam.separation(pos)
        c = r.degree<args.radius
        if np.sum(c)>0:
            print(tle.name, tle.id, np.sum(c)*args.dt)
