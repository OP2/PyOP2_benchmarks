#!/usr/bin/env python
from optparse import OptionParser
import sys
import os

meshtemplate='''
Point(1) = {0, 0, 0, %(dx)f};
Extrude {1, 0, 0} {
  Point{1}; Layers{%(layers)d};
}
Extrude {0, 1, 0} {
  Line{1}; Layers{%(layers)d};
}
'''

def generate_meshfile(name,layers):


    file(name+".geo",'w').write(meshtemplate % {'dx': 1./layers, 'layers': layers})

    os.system("gmsh -2 "+name+".geo")
    os.system("./gmsh2triangle --2d "+name+".msh")

#####################################################################
# Script starts here.

if __name__ == '__main__':
    optparser=OptionParser(usage='usage: %prog [options] <name> <layers>',
                           add_help_option=True,
                           description="""Generate the mesh files for a given"""+
                   """number of layers of elements in the channel.""")

    (options, argv) = optparser.parse_args()

    try:
        name=argv[0]
        layers=int(argv[1])
    except:
        optparser.print_help()
        sys.exit(1)

    sys.path.append(".")

    generate_meshfile(name,layers)
