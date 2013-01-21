from dolfin import *

from test_utils import *
from analytical_solution import helmholtz_initial as initial, helmholtz as analytical

parameters.form_compiler.cpp_optimize = True
parameters.krylov_solver.relative_tolerance = 1e-7

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

    M, b = assemble_system(A, RHS)
    solve(M, f.vector(), b, "cg", "jacobi")

    return sqrt(assemble((f-f_analytic)**2*dx))

if __name__ == '__main__':
    mms_main(simulation)
