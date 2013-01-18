from dolfin import *

from test_utils import *
from analytical_solution import helmholtz_initial as initial

parameters["form_compiler"]["cpp_optimize"] = True

def simulation(meshsize, degree):

    set_log_level(WARNING)

    mesh = UnitSquareMesh(meshsize, meshsize)
    mesh.init()

    # Create FunctionSpaces
    T = FunctionSpace(mesh, "CG", degree)

    # Initialise source function and previous solution function
    f = Function(T)
    f.assign(InitialCondition(initial))

    # Test and trial functions
    u, v = TrialFunction(T), TestFunction(T)

    lmbda = 1
    A = (dot(grad(v),grad(u))-lmbda*v*u)*dx
    RHS = v*f*dx

    solve(A==RHS, f, solver_parameters={"linear_solver": "cg", "preconditioner": "jacobi"})

    # Save solution in VTK format
    file = File("dolfin_helmholtz_p%d_%d.pvd" % (degree, meshsize))
    file << f

if __name__ == '__main__':
    test_main(simulation)
