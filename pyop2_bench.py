# PyOP2 vs. Fluidity vs. DOLFIN benchmark

import os
import csv
import pickle
from collections import defaultdict
from datetime import datetime
import pylab

from benchrun import Benchmark, clock
from generate_mesh import generate_meshfile

def dolfin(meshsize):
    os.system("python dolfin_adv_diff.py %d" % meshsize)

def fluidity(meshsize):
    os.system('${FLUIDITY_DIR}/bin/fluidity advection_diffusion.%d.flml' % meshsize)

def fluidity_pyop2(meshsize):
    os.system('${FLUIDITY_DIR}/bin/fluidity ufl_advection_diffusion.%d.flml' % meshsize)

def pyop2(meshsize):
    os.system('python pyop2_adv_diff.py -m mesh_%d' % meshsize)

class PyOP2Benchmark(Benchmark):
    """PyOP2 vs. Fluidity vs. DOLFIN benchmark."""

    # Execute for all combinations of these parameters
    parameters = ['version', 'meshsize']
    version = ['fluidity', 'fluidity_pyop2', 'pyop2', 'dolfin']
    meshsize = [10 * 2**i for i in range(6)]

    # Compare timings against this parameter value
    reference = ('version', 'fluidity')

    def __init__(self):
        Benchmark.__init__(self)
        self.plotdata = defaultdict(list)
        self.timestamp = datetime.today().isoformat().replace(':', '')

    def create_input(self):
        for s in self.meshsize:
            mesh = 'mesh_%d' % s
            # Generate triangle mesh
            generate_meshfile(mesh, s)
            # Generate flml
            for flml in ['advection_diffusion', 'ufl_advection_diffusion']:
                with open('%s.flml.template'%flml) as f1, open('%s.%d.flml'%(flml,s), 'w') as f2:
                    f2.write(f1.read() % {'mesh': mesh})

    def run(self, version, meshsize):
        print "Running %s with mesh size %dx%d" % (version, meshsize, meshsize)
        t1 = clock()
        globals()[version](meshsize)
        t =  clock() - t1
        elems = 2*meshsize*meshsize
        if not elems in self.plotdata['elements']:
            self.plotdata['elements'].append(elems)
        self.plotdata[version].append(t)
        return t

    def _path(self, filename):
        if not os.path.exists(self.timestamp):
            os.makedirs(self.timestamp)
        return os.path.join(self.timestamp, filename)

    def write_csv(self):
        with open(self._path('results.csv'), 'wb') as f:
            w = csv.writer(f)
            w.writerows(self.results)

    def dump(self):
        with open(self._path('results.pickle'), 'wb') as f:
            pickle.dump(self.__dict__, f)

    def load(self, filename):
        with open(filename, 'rb') as f:
            self.__dict__.update(pickle.load(f))

    def compute_speedup(self):
        for v in self.version:
            for this, base in zip(self.plotdata[v], self.plotdata['fluidity']):
                self.plotdata[v+'_speedup'].append(base/this)

    def _plot(self, fig, plot, col, legend_pos, ylabel, title, show=False):
        f = pylab.figure(fig, figsize=(8, 6), dpi=300)
        for v in self.version:
            plot(self.plotdata['elements'], self.plotdata[col(v)], label=v)
        pylab.legend(loc=legend_pos)
        pylab.xlabel('Number of elements in the mesh')
        pylab.ylabel(ylabel)
        pylab.title(title)
        pylab.grid(True)
        pylab.savefig(self._path('%s.svg' % fig), orientation='landscape', format='svg', transparent=True)
        if show:
            pylab.show()
        pylab.close(f)

    def plot_result(self):
        title = 'Benchmark of an advection-diffusion problem for 100 time steps'
        for fig, pl in zip(('linear', 'semilogx'), (pylab.plot, pylab.semilogx)):
            self._plot('runtime_'+fig, pl, lambda x: x, 'upper left', 'Overall runtime in seconds', title)
            self._plot('speedup_'+fig, pl, lambda x: x+'_speedup', 'lower right', 'Relative speedup over Fluidity baseline', title)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description=PyOP2Benchmark.__doc__)
    parser.add_argument('-r', '--run', action='store_true', help='run benchmarks')
    parser.add_argument('-p', '--plot', action='store_true', help='plot results')
    parser.add_argument('-d', '--dump', action='store_true', help='Pickle dump')
    parser.add_argument('-c', '--csv', action='store_true', help='Dump to CSV')
    parser.add_argument('-l', '--load', help='Pickle load from file')
    args = parser.parse_args()

    b = PyOP2Benchmark()
    if args.load and os.path.exists(args.load):
        b.load(args.load)
    if args.run:
        b.create_input()
        b.print_result()
        b.compute_speedup()
    if args.dump:
        b.dump()
    if args.csv:
        b.write_csv()
    if args.plot:
        b.plot_result()
