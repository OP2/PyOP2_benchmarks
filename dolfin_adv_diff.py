import sys

from dolfin import *
from time import time
from analytical_solution import advection_diffusion as val

parameters["form_compiler"]["cpp_optimize"] = True

class InitialCondition(Expression):
    def __init__(self, fn):
        self._fn = fn

    def eval(self, values, x):
        values[0] = self._fn(x, 0.1)

def simulation(D, t, dt, endtime, mesh, initial):

    set_log_level(ERROR)

    # Added due to mesh not conforming to UFC numbering, allegedly
    mesh.order()

    # Create FunctionSpaces
    T = FunctionSpace(mesh, "CG", 1)

    # Create velocity Function
    velocity = Constant( (1.0, 0.0) )

    concentration = InitialCondition(initial)

    # Initialise source function and previous solution function
    u0 = Function(T)
    u1 = Function(T)
    u1.assign(concentration)

    # Test and trial functions
    u, v = TrialFunction(T), TestFunction(T)

    # Advection
    Mass = v*u*dx
    adv_rhs = (v*u0+dt*dot(grad(v),velocity)*u0)*dx

    #Diffusion
    d = -dt*D*dot(grad(v),grad(u))*dx
    diff_matrix = Mass - 0.5*d
    diff_rhs = action(Mass + 0.5*d, u0)

    M_a, b_a, M_d, b_d = None, None, None, None
    # Time-stepping
    first = True
    while t < endtime:

        if first:
            t1 = time()
            first = False
        # Copy soln from prev.
        u0.assign(u1)

        # Advection
        M_a = assemble(Mass, tensor=M_a, reset_sparsity=(M_a is None))
        b_a = assemble(adv_rhs, tensor=b_a, reset_sparsity=(b_a is None))
        solve(M_a, u1.vector(), b_a, "cg", "jacobi")

        # Copy solution from advection
        u0.assign(u1)

        # Diffusion
        M_d = assemble(diff_matrix, tensor=M_d, reset_sparsity=(M_d is None))
        b_d = assemble(diff_rhs, tensor=b_d, reset_sparsity=(b_d is None))
        solve(M_d, u1.vector(), b_d, "cg", "jacobi")

        # Next timestep
        t += dt
    
    if MPI.process_number() == 0:
        print "/fluidity :: %f" % (time()-t1)

def run(meshsize):
    from parameters import diffusivity, current_time, dt, endtime
    mesh = Mesh('meshes/square.%s.xml' % meshsize)
    mesh.init()

    simulation(diffusivity, current_time, dt, endtime, mesh, val)

if __name__ == '__main__':
    run(sys.argv[1])
