"""Concatenate cProfile results from multiple files.

Usage: run PATTERN FILE [mpi]
"""

import os
import sys
from glob import glob
from pstats import Stats

def concat(pattern, outfile, mpi=None):
    if mpi:
        from mpi4py import MPI
        pattern = pattern % MPI.COMM_WORLD.rank
        outfile = outfile % MPI.COMM_WORLD.rank
    files = glob(pattern)
    if files:
        s = Stats(files[0])
        for f in files[1:]: s.add(f)
        s.dump_stats(outfile)
        for f in files:
            os.remove(f)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(1)
    concat(sys.argv[1], sys.argv[2], len(sys.argv) == 4)
