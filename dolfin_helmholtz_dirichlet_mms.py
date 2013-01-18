from dolfin import *

from test_utils import *
from analytical_solution import helmholtz_initial_sin as initial, helmholtz_sin as analytical

parameters["form_compiler"]["cpp_optimize"] = True

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

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(T, u0, lambda x:
            x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
            x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)

    solve(A==RHS, f, bcs=bc, solver_parameters={"linear_solver": "cg", "preconditioner": "jacobi"})

    return sqrt(assemble((f-f_analytic)**2*dx))

if __name__ == '__main__':
    mms_main(simulation)
