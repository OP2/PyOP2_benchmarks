#!/usr/bin/env python

"Generate the mesh files for a given number of layers of elements in the channel."

from argparse import ArgumentParser
import sys
import subprocess

meshtemplate='''
Point(1) = {0, 0, 0, %(dx)f};
Extrude {1, 0, 0} {
  Point{1}; Layers{%(layers)d};
}
Extrude {0, 1, 0} {
  Line{1}; Layers{%(layers)d};
}
'''

def generate_meshfile(name, layers, capture=False):

    with open(name+".geo",'w') as f:
        f.write(meshtemplate % {'dx': 1./layers, 'layers': layers})

    cmd = "gmsh -2 "+name+".geo && ./gmsh2triangle --2d "+name+".msh"
    if capture:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    return subprocess.call(cmd, shell=True)

#####################################################################
# Script starts here.

if __name__ == '__main__':
    parser=ArgumentParser(add_help=True, description=__doc__)
    parser.add_argument('name', help='Base name for the generated output files')
    parser.add_argument('layers', type=int, help='Number of layers to generate')
    args = parser.parse_args()

    sys.exit(generate_meshfile(args.name, args.layers))
