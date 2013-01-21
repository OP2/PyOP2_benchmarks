from dolfin import *

from test_utils import *
from analytical_solution import helmholtz_initial_sin as initial

parameters.form_compiler.cpp_optimize = True
parameters.krylov_solver.relative_tolerance = 1e-7

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

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(T, u0, lambda x:
            x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
            x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)

    M, b = assemble_system(A, RHS, bc)
    solve(M, f.vector(), b, "cg", "jacobi")

    # Save solution in VTK format
    file = File("dolfin_helmholtz_dirichlet_p%d_%d.pvd" % (degree, meshsize))
    file << f

if __name__ == '__main__':
    test_main(simulation)
