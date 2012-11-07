import logging
import os
import pylab

from benchrun import clock
from pyop2_bench import PyOP2Benchmark
from generate_mesh import generate_meshfile
import parameters

class AdvDiffBenchmark(PyOP2Benchmark):
    """PyOP2 vs. Fluidity vs. DOLFIN benchmark."""

    # Execute for all combinations of these parameters
    parameters = ['version', 'meshsize']
    version = ['fluidity', 'fluidity_pyop2', 'pyop2', 'dolfin']
    meshsize = [10 * 2**i for i in range(6)]

    # Compare timings against this parameter value
    reference = ('version', 'fluidity')

    def dolfin(self, meshsize):
        self.logged_call("python dolfin_adv_diff.py %d" % meshsize)

    def fluidity(self, meshsize):
        self.logged_call('${FLUIDITY_DIR}/bin/fluidity advection_diffusion.%d.flml' % meshsize)

    def fluidity_pyop2(self, meshsize):
        self.logged_call('${FLUIDITY_DIR}/bin/fluidity ufl_advection_diffusion.%d.flml' % meshsize)

    def pyop2(self, meshsize):
        self.logged_call('python pyop2_adv_diff.py -m mesh_%d' % meshsize)

    def create_input(self):
        for s in self.meshsize:
            mesh = 'mesh_%d' % s
            # Generate triangle mesh
            self.log(generate_meshfile(mesh, s, capture=True))
            # Generate flml
            for flml in ['advection_diffusion', 'ufl_advection_diffusion']:
                with open('%s.flml.template'%flml) as f1, open('%s.%d.flml'%(flml,s), 'w') as f2:
                    f2.write(f1.read() % {
                        'mesh': mesh,
                        'diffusivity': parameters.diffusivity,
                        'current_time': parameters.current_time,
                        'dt': parameters.dt,
                        'endtime': (parameters.endtime - parameters.dt)
                        })

    def run(self, version, meshsize):
        self.log("Running %s with mesh size %dx%d" % (version, meshsize, meshsize))
        t1 = clock()
        self.__getattribute__(version)(meshsize)
        t =  clock() - t1
        elems = 2*meshsize*meshsize
        if not elems in self.plotdata['elements']:
            self.plotdata['elements'].append(elems)
        self.plotdata[version].append(t)
        return t

    def time_all(self):
        self.log('PyOP2 version:')
        self.logged_call('GIT_DIR=${PYOP2_DIR}/.git git rev-parse HEAD')
        self.log('Fluidity version:')
        self.logged_call('bzr revision-info -d ${FLUIDITY_DIR}')
        self.log('UFL version:')
        self.logged_call('bzr revision-info -d ${UFL_DIR}')
        self.log('FFC version:')
        self.logged_call('bzr revision-info -d ${FFC_DIR}')
        self.log('DOLFIN version:')
        self.logged_call('bzr revision-info -d ${DOLFIN_DIR}')
        super(AdvDiffBenchmark, self).time_all()

    def plot_result(self):
        title = 'Benchmark of an advection-diffusion problem for 100 time steps'
        for fig, pl in zip(('linear', 'semilogx'), (pylab.plot, pylab.semilogx)):
            self._plot('runtime_'+fig, pl, lambda x: x, 'upper left', 'Overall runtime in seconds', title)
            self._plot('speedup_'+fig, pl, lambda x: x+'_speedup', 'lower right', 'Relative speedup over Fluidity baseline', title)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description=AdvDiffBenchmark.__doc__)
    parser.add_argument('-r', '--run', action='store_true', help='run benchmarks')
    parser.add_argument('-p', '--plot', action='store_true', help='plot results')
    parser.add_argument('-d', '--dump', action='store_true', help='Pickle dump')
    parser.add_argument('-c', '--csv', action='store_true', help='Dump to CSV')
    parser.add_argument('-l', '--load', help='Pickle load from file')
    parser.add_argument('-q', '--quiet', help='Only print errors and warnings')
    args = parser.parse_args()

    b = AdvDiffBenchmark()
    logging.getLogger().setLevel(logging.WARN if args.quiet else logging.INFO)
    if args.load and os.path.exists(args.load):
        b.load(args.load)
    if args.run:
        b.create_input()
        b.time_all()
        b.sort_results()
        b.print_result()
        b.compute_speedup()
    if args.dump:
        b.dump()
    if args.csv:
        b.write_csv()
    if args.plot:
        b.plot_result()
