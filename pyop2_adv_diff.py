"""
This demo solves an advection-diffusion equation on a domain read in from a
triangle file. It requires the pyop2 branch of ffc, which can be obtained
with:

bzr branch lp:~mapdes/ffc/pyop2

This may also depend on development trunk versions of other FEniCS programs.
"""

import numpy as np

from benchrun import clock
from pyop2 import op2, utils
from pyop2.ffc_interface import compile_form
from triangle_reader import read_triangle
from ufl import *

def run(diffusivity, current_time, dt, endtime, **kwargs):
    op2.init(**kwargs)

    # Set up finite element problem

    T = FiniteElement("Lagrange", "triangle", 1)
    V = VectorElement("Lagrange", "triangle", 1)

    p=TrialFunction(T)
    q=TestFunction(T)
    t=Coefficient(T)
    u=Coefficient(V)

    diffusivity = 0.1

    M=p*q*dx

    adv_rhs = (q*t+dt*dot(grad(q),u)*t)*dx

    d=-dt*diffusivity*dot(grad(q),grad(p))*dx

    diff_matrix=M-0.5*d
    diff_rhs=action(M+0.5*d,t)

    # Generate code for mass and rhs assembly.

    mass        = compile_form(M,           "mass")[0]
    adv_rhs     = compile_form(adv_rhs,     "adv_rhs")[0]
    diff_matrix = compile_form(diff_matrix, "diff_matrix")[0]
    diff_rhs    = compile_form(diff_rhs,    "diff_rhs")[0]

    # Set up simulation data structures

    valuetype=np.float64

    nodes, coords, elements, elem_node = read_triangle(kwargs['mesh'])
    num_nodes = nodes.size

    sparsity = op2.Sparsity((elem_node, elem_node), 1, "sparsity")
    mat = op2.Mat(sparsity, valuetype, "mat")

    tracer_vals = np.asarray([0.0]*num_nodes, dtype=valuetype)
    tracer = op2.Dat(nodes, 1, tracer_vals, valuetype, "tracer")

    b_vals = np.asarray([0.0]*num_nodes, dtype=valuetype)
    b = op2.Dat(nodes, 1, b_vals, valuetype, "b")

    velocity_vals = np.asarray([1.0, 0.0]*num_nodes, dtype=valuetype)
    velocity = op2.Dat(nodes, 2, velocity_vals, valuetype, "velocity")

    # Set initial condition

    i_cond_code="""
    void i_cond(double *c, double *t)
    {
      double i_t = 0.01; // Initial time
      double A   = 0.1; // Normalisation
      double D   = 0.1; // Diffusivity
      double pi  = 3.141459265358979;
      double x   = c[0]-0.5;
      double y   = c[1]-0.5;
      double r   = sqrt(x*x+y*y);

      if (r<0.25)
        *t = A*(exp((-(r*r))/(4*D*i_t))/(4*pi*D*i_t));
      else
        *t = 0.0;
    }
    """

    i_cond = op2.Kernel(i_cond_code, "i_cond")

    op2.par_loop(i_cond, nodes,
                 coords(op2.IdentityMap, op2.READ),
                 tracer(op2.IdentityMap, op2.WRITE))

    zero_dat_code="""
    void zero_dat(double *dat)
    {
      *dat = 0.0;
    }
    """

    zero_dat = op2.Kernel(zero_dat_code, "zero_dat")

    # Assemble and solve

    have_advection = True
    have_diffusion = True

    t1 = clock()

    while current_time < endtime:

        # Advection

        if have_advection:
            mat.zero()

            op2.par_loop(mass, elements(3,3),
                         mat((elem_node[op2.i[0]], elem_node[op2.i[1]]), op2.INC),
                         coords(elem_node, op2.READ))

            op2.par_loop(zero_dat, nodes,
                         b(op2.IdentityMap, op2.WRITE))

            op2.par_loop(adv_rhs, elements(3),
                         b(elem_node[op2.i[0]], op2.INC),
                         coords(elem_node, op2.READ),
                         tracer(elem_node, op2.READ),
                         velocity(elem_node, op2.READ))

            op2.solve(mat, tracer, b)

        # Diffusion

        if have_diffusion:
            mat.zero()

            op2.par_loop(diff_matrix, elements(3,3),
                         mat((elem_node[op2.i[0]], elem_node[op2.i[1]]), op2.INC),
                         coords(elem_node, op2.READ))

            op2.par_loop(zero_dat, nodes,
                         b(op2.IdentityMap, op2.WRITE))

            op2.par_loop(diff_rhs, elements(3),
                         b(elem_node[op2.i[0]], op2.INC),
                         coords(elem_node, op2.READ),
                         tracer(elem_node, op2.READ))

            op2.solve(mat, tracer, b)


        current_time += dt

    runtime = clock() - t1
    print "/fluidity :: %f" % runtime

if __name__ == '__main__':
    from parameters import *
    parser = utils.parser(group=True, description="PyOP2 P1 advection-diffusion demo")
    parser.add_argument('-m', '--mesh',
                        help='Base name of triangle mesh (excluding the .ele or .node extension)')
    opt = vars(parser.parse_args())
    run(diffusivity, current_time, dt, endtime, **opt)
