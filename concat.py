"""Usage: run PATTERN FILE"""

import sys
from glob import glob
from pstats import Stats

if len(sys.argv) != 3:
    print __doc__
    sys.exit(1)
files = glob(sys.argv[1])
s = Stats(files[0])
for f in files[1:]: s.add(f)
s.dump_stats(sys.argv[2])
