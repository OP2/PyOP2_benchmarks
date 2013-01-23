from math import sqrt
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
    version = ['fluidity', 'fluidity_pyop2_seq', 'pyop2_seq']
    meshsize = [int(10 * sqrt(2)**i) for i in range(11)]

    # Compare timings against this parameter value
    reference = ('version', 'fluidity')

    plotstyle = dict(zip(version, ['k-o', 'g-s', 'r-d', 'b-^']))

    def dolfin(self, meshsize):
        self.logged_call(self.mpicmd+"python dolfin_adv_diff.py %d" % meshsize)

    def fluidity(self, meshsize):
        self.logged_call(self.mpicmd+'${FLUIDITY_DIR}/bin/fluidity flmls/advection_diffusion.%d.flml' % meshsize)

    def fluidity_pyop2_seq(self, meshsize):
        self.logged_call('${FLUIDITY_DIR}/bin/fluidity flmls/ufl_advection_diffusion.sequential.%d.flml' % meshsize)

    def fluidity_pyop2_openmp(self, meshsize):
        self.logged_call('${FLUIDITY_DIR}/bin/fluidity flmls/ufl_advection_diffusion.openmp.%d.flml' % meshsize)

    def fluidity_pyop2_cuda(self, meshsize):
        self.logged_call('${FLUIDITY_DIR}/bin/fluidity flmls/ufl_advection_diffusion.cuda.%d.flml' % meshsize)

    def pyop2_seq(self, meshsize):
        self.logged_call('python pyop2_adv_diff.py -m meshes/mesh_%d -b sequential' % meshsize)

    def pyop2_openmp(self, meshsize):
        self.logged_call('python pyop2_adv_diff.py -m meshes/mesh_%d -b openmp' % meshsize)

    def pyop2_cuda(self, meshsize):
        self.logged_call('python pyop2_adv_diff.py -m meshes/mesh_%d -b cuda' % meshsize)

    def __init__(self, np=1, message=''):
        super(AdvDiffBenchmark, self).__init__()
        self.log("Running versions: %s" % AdvDiffBenchmark.version)
        self.log("Reference %s: %s" % AdvDiffBenchmark.reference)
        self.np=np
        self.message=message
        self.mpicmd = 'mpirun --bycore --bysocket --bind-to-socket --bind-to-core -np %d ' % np if np > 1 else ''
        self.plotlabels = {
                'fluidity': 'Fluidity (cores: %d)' % np,
                'fluidity_pyop2_seq': 'Fluidity-PyOP2 (backend: sequential)',
                'pyop2_seq': 'PyOP2 (backend: sequential)',
                'dolfin': 'DOLFIN (cores: %d)' % np
                }

    def create_input(self):
        if not os.path.exists('meshes'):
            os.makedirs('meshes')
        if not os.path.exists('flmls'):
            os.makedirs('flmls')
        for s in self.meshsize:
            mesh = os.path.join('meshes','mesh_%d' % s)
            # Generate triangle mesh
            self.log(generate_meshfile(mesh, s, capture=True))
            # Decompose the mesh if running in parallel
            if self.np > 1:
                self.logged_call('${FLUIDITY_DIR}/bin/fldecomp -m triangle -n %d %s' \
                        % (self.np, mesh))
            # Generate flml
            for backend in ['sequential', 'openmp', 'cuda']:
                with open('ufl_advection_diffusion.flml.template') as f1, \
                        open(os.path.join('flmls','ufl_advection_diffusion.%s.%d.flml'%(backend,s)), 'w') as f2:
                    f2.write(f1.read() % {
                        'mesh': mesh,
                        'diffusivity': parameters.diffusivity,
                        'current_time': parameters.current_time,
                        'dt': parameters.dt,
                        'endtime': (parameters.endtime - parameters.dt),
                        'backend': backend
                        })
            with open('advection_diffusion.flml.template') as f1, \
                        open(os.path.join('flmls','advection_diffusion.%d.flml'%s), 'w') as f2:
                    f2.write(f1.read() % {
                        'mesh': mesh,
                        'diffusivity': parameters.diffusivity,
                        'current_time': parameters.current_time,
                        'dt': parameters.dt,
                        'endtime': (parameters.endtime - parameters.dt),
                        'backend': backend
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
        if self.message:
            self.log(self.message+'\n')
        self.log('PyOP2 version:')
        self.logged_call('git --git-dir=%s/.git rev-parse HEAD' % os.environ['PYOP2_DIR'])
        self.log('Fluidity version:')
        self.logged_call('bzr revision-info -d %s' % os.environ['FLUIDITY_DIR'])
        self.log('UFL version:')
        self.logged_call('bzr revision-info -d %s' % os.environ['UFL_DIR'])
        self.log('FFC version:')
        self.logged_call('bzr revision-info -d %s' % os.environ['FFC_DIR'])
        #self.log('DOLFIN version:')
        #self.logged_call('bzr revision-info -d %s' % os.environ['DOLFIN_DIR'])
        super(AdvDiffBenchmark, self).time_all()

    def plot_result(self, format='svg'):
        title = 'Benchmark of an advection-diffusion problem for 100 time steps'
        title = self.message or 'Benchmark of an advection-diffusion problem for 100 time steps'
        for fig, pl in zip(('linear', 'semilogx'), (pylab.plot, pylab.semilogx)):
            self._plot('runtime_'+fig, pl, lambda x: x, 'upper left', 'Overall runtime in seconds', title, format)
            self._plot('speedup_'+fig, pl, lambda x: x+'_speedup', 'lower right', 'Relative speedup over Fluidity baseline', title, format)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description=AdvDiffBenchmark.__doc__)
    parser.add_argument('-r', '--run', action='store_true', help='run benchmarks')
    parser.add_argument('-p', '--plot', action='store_true', help='plot results')
    parser.add_argument('-d', '--dump', action='store_true', help='Pickle dump')
    parser.add_argument('-c', '--csv', action='store_true', help='Dump to CSV')
    parser.add_argument('-v', '--versions', nargs='+', help='Versions to run')
    parser.add_argument('-a', '--reference', help='Version to use as reference')
    parser.add_argument('-l', '--load', help='Pickle load from file')
    parser.add_argument('-q', '--quiet', help='Only print errors and warnings')
    parser.add_argument('-s', '--create-input', action='store_true',
            help='Do not generate input files')
    parser.add_argument('-n', '-np', type=int, default=1,
            help='Number of MPI process (Fluidity, DOLFIN)')
    parser.add_argument('-m', '--message', default='',
            help='Message, added to the log output')
    args = parser.parse_args()

    if args.versions:
        AdvDiffBenchmark.version = args.versions
    if args.reference:
        AdvDiffBenchmark.reference = ('version', args.reference)

    if args.plot:
        print "going to plot."

    logging.getLogger().setLevel(logging.WARN if args.quiet else logging.INFO)
    b = AdvDiffBenchmark(args.n, args.message)
    if args.load and os.path.exists(args.load):
        b.load(args.load)
    if args.create_input:
        b.create_input()
    if args.run:
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
