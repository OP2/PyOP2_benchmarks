from dolfin import Expression
import numpy as np

class InitialCondition(Expression):
    def __init__(self, fn, n=8, t=0.0):
        self._fn = fn
        self._n = n
        self._t = t

    def eval(self, values, x):
        values[0] = self._fn(x, self._t, self._n)

def test_main(f):
    import sys
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    if len(sys.argv) > 1:
        f(int(sys.argv[1]), d)
    else:
        for d in range(1,3):
            for i in range(4):
                f(30*2**i, d)

def mms_main(f):
    for d in range(1,3):
        err = np.array([f(30*2**i, d) for i in range(4)])
        conv = np.log2(err[:-1]/err[1:])
        print "P%d: L2 error norms: %s, convergence order: %s" % (d, err, conv)
