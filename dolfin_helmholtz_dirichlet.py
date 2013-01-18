from dolfin import *
from analytical_solution import helmholtz_initial_sin as initial

parameters["form_compiler"]["cpp_optimize"] = True

class InitialCondition(Expression):
    def __init__(self, fn):
        self._fn = fn

    def eval(self, values, x):
        values[0] = self._fn(x, 0.0)

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

    solve(A==RHS, f, bcs=bc, solver_parameters={"linear_solver": "cg", "preconditioner": "jacobi"})

    # Save solution in VTK format
    file = File("dolfin_helmholtz_dirichlet_p%d_%d.pvd" % (degree, meshsize))
    file << f

if __name__ == '__main__':
    import sys
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    if len(sys.argv) > 1:
        simulation(int(sys.argv[1]), d)
    else:
        for d in range(1,3):
            for i in range(4):
                simulation(30*2**i, d)
