"""PyOP2 sparsity building benchmark"""

import numpy as np
from time import clock

from pyop2 import op2, utils

def time_sparsity(n):
    nodes = op2.Set(n)
    cells = op2.Set(n-1)
    m = op2.Map(cells, nodes, 2, np.concatenate((np.arange(n-1), np.arange(1,n))))
    t = clock()
    s = op2.Sparsity((m,m), 1)
    return clock() - t

if __name__=='__main__':
    import pylab
    parser = utils.parser(group=True, description=__doc__)
    op2.init(**vars(parser.parse_args()))
    n = [2**i for i in range(10,18)]
    t = [time_sparsity(i) for i in n]
    pylab.plot(n,t)
    pylab.xlabel('Sparsity size')
    pylab.ylabel('Sparsity build time (s)')
    pylab.show()
