#!/usr/bin/env python

"Generate the mesh files for a given number of layers of elements in the channel."

from argparse import ArgumentParser
import os
import sys
import subprocess


def generate_trianglefile(mesh, size, capture=False, reorder=True, move=True, cwd=None):
    if capture:
        runcmd = lambda cmd: subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, cwd=cwd)
    else:
        runcmd = lambda cmd: subprocess.call(cmd, shell=True, cwd=cwd)
    generate = runcmd("triangle -e -a%s %s" % (size, 'square.poly'))

    if reorder:
        fluidity = runcmd("${HILBERT_DIR}/bin/fluidity " + 'reorder_mesh.flml')
        generate = '\n'.join([generate, fluidity])
        if move:
            movele = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.ele %s.ele" % mesh)
            movedge = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.edge %s.edge" % mesh)
            movnode = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.node %s.node" % mesh)
            generate = '\n'.join([generate, movele, movedge, movnode])
    elif move:
        movele = runcmd("mv square.1.ele %s.ele" % mesh)
        movedge = runcmd("mv square.1.edge %s.edge" % mesh)
        movnode = runcmd("mv square.1.node %s.node" % mesh)
        generate = '\n'.join([generate, movele, movedge, movnode])
    runcmd("rm -f square.1.*")
    if capture:
        return generate

#####################################################################
# Script starts here.

if __name__ == '__main__':
    parser=ArgumentParser(add_help=True, description=__doc__)
    parser.add_argument('size', type=int, help='Size to pass to triangle')
    args = parser.parse_args()

    sys.exit(generate_meshfile(args.name, args.layers))
