from dolfin import *
import numpy as np

from analytical_solution import helmholtz_initial as initial, helmholtz as analytical

parameters["form_compiler"]["cpp_optimize"] = True

class InitialCondition(Expression):
    def __init__(self, fn, n=8):
        self._fn = fn
        self._n = n

    def eval(self, values, x):
        values[0] = self._fn(x, 0.0, self._n)

def simulation(meshsize, degree, n=8):

    set_log_level(WARNING)

    mesh = UnitSquareMesh(meshsize, meshsize)
    mesh.init()

    # Create FunctionSpaces
    T = FunctionSpace(mesh, "CG", degree)

    # Initialise source function and previous solution function
    f = Function(T)
    f.assign(InitialCondition(initial, n))
    f_analytic = Function(T)
    f_analytic.assign(InitialCondition(analytical, n))

    # Test and trial functions
    u, v = TrialFunction(T), TestFunction(T)

    lmbda = 1
    A = (dot(grad(v),grad(u))-lmbda*v*u)*dx
    RHS = v*f*dx

    solve(A==RHS, f, solver_parameters={"linear_solver": "cg", "preconditioner": "jacobi"})

    return sqrt(assemble((f-f_analytic)**2*dx))

if __name__ == '__main__':
    for d in range(1,3):
        err = np.array([simulation(30*2**i, d) for i in range(4)])
        conv = np.log2(err[:-1]/err[1:])
        print "P%d: L2 error norms: %s, convergence order: %s" % (d, err, conv)
