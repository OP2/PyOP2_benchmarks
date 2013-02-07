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

    def dolfin(self, meshsize):
        return self.logged_call_with_time(self.mpicmd+"python dolfin_adv_diff.py %d" % meshsize)

    def fluidity(self, meshsize):
        return self.logged_call_with_time('${FLUIDITY_DIR}/bin/fluidity -p flmls/advection_diffusion.%d.flml' % meshsize)

    def fluidity_mpi(self, meshsize):
        return self.logged_call_with_time('OMP_NUM_THREADS=1 '+self.mpicmd+'${FLUIDITY_DIR}/bin/fluidity -p flmls/advection_diffusion.%d.flml' % meshsize)

    def fluidity_pyop2_seq(self, meshsize):
        return self.flufl_call_with_time('${FLUIDITY_DIR}/bin/fluidity -p flmls/ufl_advection_diffusion.sequential.%d.flml' % meshsize)

    def fluidity_pyop2_openmp(self, meshsize):
        return self.flufl_call_with_time('${FLUIDITY_DIR}/bin/fluidity -p flmls/ufl_advection_diffusion.openmp.%d.flml' % meshsize)

    def fluidity_pyop2_cuda(self, meshsize):
        return self.flufl_call_with_time('${FLUIDITY_DIR}/bin/fluidity -p flmls/ufl_advection_diffusion.cuda.%d.flml' % meshsize)

    def pyop2_seq(self, meshsize):
        return self.logged_call_with_time('python pyop2_adv_diff.py -m meshes/mesh_%d -b sequential' % meshsize)

    def pyop2_openmp(self, meshsize):
        return self.logged_call_with_time('python pyop2_adv_diff.py -m meshes/mesh_%d -b openmp' % meshsize)

    def pyop2_cuda(self, meshsize):
        return self.logged_call_with_time('python pyop2_adv_diff.py -m meshes/mesh_%d -b cuda' % meshsize)

    def __init__(self, np=1, message='', version=None, reference=None, meshsize=None, parameters=None):
        # Execute for all combinations of these parameters
        #self.meshsize = meshsize or [101, 142, 174, 201, 225, 246, 266, 284, 301, 317, 333, 347, 362, 375, 388]
        #self.meshsize = meshsize or [101, 174, 225, 266, 301, 333, 362, 388]
        self.meshsize = [317, 448, 549, 633]#, 708]
        #[int(10 * sqrt(2)**i) for i in range(11)]
        parameters = parameters or ['version', 'meshsize']
        self.version = version or ['fluidity', 'fluidity_pyop2_seq', 'pyop2_seq']
        self.primed = []
        # Compare timings against this parameter value
        self.reference = ('version', reference or 'fluidity')
        super(AdvDiffBenchmark, self).__init__(parameters)

        self.log("Running versions: %s" % self.version)
        self.log("Reference %s: %s" % self.reference)

        self.np=np
        self.message=message
        self.mpicmd = 'mpiexec ' if np>1 else ''

        self.plotstyle = dict(zip(self.version, ['k-o', 'g-s', 'r-d', 'b-^']))
        self.plotlabels = {
                'fluidity': 'Fluidity (cores: 1)',
                'fluidity_mpi': 'Fluidity (cores: %d)' % self.np,
                'fluidity_pyop2_seq': 'Fluidity-PyOP2 (backend: sequential)',
                'fluidity_pyop2_cuda': 'Fluidity-PyOP2 (backend: cuda)',
                'fluidity_pyop2_openmp': 'Fluidity-PyOP2 (backend: openmp)',
                'pyop2_seq': 'PyOP2 (backend: sequential)',
                'pyop2_openmp': 'PyOP2 (backend: OpenMP)',
                'pyop2_cuda': 'PyOP2 (backend: CUDA)',
                'dolfin': 'DOLFIN (cores: %d)' % self.np
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


    def dry_run(self, version, meshsize):
        if version not in self.primed:
            self.log("Performing dry run of %s using mesh size %dx%d" % (version, meshsize, meshsize))
            self.__getattribute__(version)(meshsize)
            self.primed.append(version)

    def run(self, version, meshsize):
        self.dry_run(version, meshsize)
        self.log("Running %s with mesh size %dx%d" % (version, meshsize, meshsize))
        t = self.__getattribute__(version)(meshsize)
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

    logging.getLogger().setLevel(logging.WARN if args.quiet else logging.INFO)
    b = AdvDiffBenchmark(args.n, args.message, args.versions, args.reference)
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
