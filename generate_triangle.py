#!/usr/bin/env python

"Generate the mesh files for a given number of layers of elements in the channel."

from argparse import ArgumentParser
import sys
import subprocess

def runcmd(cmd, capture):
    if capture:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    return subprocess.call(cmd, shell=True)

def generate_trianglefile(size, capture=False, reorder=True, move=True):
    generate = runcmd("triangle -e -a"+size+" square.poly", capture)
    if reorder:
        fluidity = runcmd("${HILBERT_DIR}/bin/fluidity reorder_mesh.flml", capture)
        generate = '\n'.join([generate, fluidity])
        if move:
            movele = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.ele meshes/square."+size+".ele", capture)
            movedge = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.edge meshes/square."+size+".edge", capture)
            movnode = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.node meshes/square."+size+".node", capture)
            #movpoly = runcmd("mv reorder_mesh_CoordinateMesh_0_checkpoint.poly meshes/square."+size+".poly", capture)
            generate = '\n'.join([generate, movele, movedge, movnode])
    elif move:
        movele = runcmd("mv square.1.ele meshes/square."+size+".ele", capture)
        movedge = runcmd("mv square.1.edge meshes/square."+size+".edge", capture)
        movnode = runcmd("mv square.1.node meshes/square."+size+".node", capture)
        #movpoly = runcmd("mv square.1.poly meshes/square."+size+".poly", capture)
        generate = '\n'.join([generate, movele, movedge, movnode])
    if capture:
        return generate

#####################################################################
# Script starts here.

if __name__ == '__main__':
    parser=ArgumentParser(add_help=True, description=__doc__)
    parser.add_argument('size', type=int, help='Size to pass to triangle')
    args = parser.parse_args()

    sys.exit(generate_meshfile(args.name, args.layers))
