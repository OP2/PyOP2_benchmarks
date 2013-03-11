"""Usage: run PATTERN FILE"""

import os
import sys
from glob import glob
from pstats import Stats

def concat(pattern, outfile):
    files = glob(pattern)
    s = Stats(files[0])
    for f in files[1:]: s.add(f)
    s.dump_stats(outfile)
    os.remove(files)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print __doc__
        sys.exit(1)
    concat(sys.argv[1], sys.argv[2])
