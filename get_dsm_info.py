#! /usr/bin/env python

# Command line tool to call on a pairname dir
import os
import sys

from osgeo import gdal, ogr, osr

import dsm_info

def main(argv=None):

    if len(sys.argv[1:]) == 1 and os.path.exists(sys.argv[1]):
        imageDir = sys.argv[1]

    else:
        sys.exit("Usage: %s [pairname_dir , i.e., /full/path/<SENSORID>_YYYYMMDD_<catid1>_<catid2> ]" % os.path.basename(sys.argv[0]))

    print(dsm_info.main(imageDir))

if __name__ == '__main__':
    main()
