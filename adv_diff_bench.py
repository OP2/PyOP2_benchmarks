from glob import glob
from math import sqrt
import logging
import os
import pylab

from benchrun import clock
from pyop2_bench import PyOP2Benchmark
from generate_mesh import generate_meshfile
from generate_triangle import generate_trianglefile
import parameters

class AdvDiffBenchmark(PyOP2Benchmark):
    """PyOP2 vs. Fluidity vs. DOLFIN benchmark."""

    def dolfin(self, meshsize):
        """DOLFIN advection-diffusion benchmark (MPI parallel)"""
        cmd = [self.mpicmd, "python dolfin_adv_diff.py %s" % meshsize]
        return self.logged_call_with_time(cmd)

    def fluidity(self, meshsize):
        """Fluidity advection-diffusion benchmark (sequential)"""
        cmd = [self.flcmd, '-p %s/advection_diffusion.%s.flml' % (self.flmldir, meshsize)]
        return self.logged_call_with_time(cmd)

    def fluidity_mpi(self, meshsize):
        """Fluidity advection-diffusion benchmark (MPI parallel)"""
        cmd = [self.mpicmd, self.flcmd, '-p %s/advection_diffusion.%s.flml'
                % (self.flmldir, meshsize)]
        return self.logged_call_with_time(cmd)

    def _fluidity_pyop2_call(self, backend, meshsize, mpi=''):
        cmd = [mpi, self.flcmd, '-p %s/ufl_advection_diffusion%s.%s.%s.flml'
                % (self.flmldir, self.profile, backend, meshsize)]
        time = self.flufl_call_with_time(cmd)
        if self.profile:
            pattern = 'ufl_advection_diffusion.%s.*.%%d.cprofile.part' % backend
            outfile = self._path('ufl_advection_diffusion.%s.%s.%%d.cprofile' % (backend, meshsize))
            self.logged_call([mpi, 'python concat.py', "'"+pattern+"'", outfile, 'mpi'])
        return time

    def fluidity_pyop2_seq(self, meshsize):
        """Fluidity-PyOP2 advection-diffusion benchmark (sequential backend)"""
        return self._fluidity_pyop2_call('sequential', meshsize)

    def fluidity_pyop2_openmp(self, meshsize):
        """Fluidity-PyOP2 advection-diffusion benchmark (OpenMP backend)"""
        return self._fluidity_pyop2_call('openmp', meshsize)

    def fluidity_pyop2_cuda(self, meshsize):
        """Fluidity-PyOP2 advection-diffusion benchmark (CUDA backend)"""
        return self._fluidity_pyop2_call('cuda', meshsize)

    def fluidity_pyop2_mpi(self, meshsize):
        """Fluidity-PyOP2 advection-diffusion benchmark (MPI parallel)"""
        return self._fluidity_pyop2_call('sequential', meshsize, self.mpicmd)

    def fluidity_pyop2_mpi_openmp(self, meshsize):
        """Fluidity-PyOP2 advection-diffusion benchmark (MPI + OpenMP parallel)"""
        return self._fluidity_pyop2_call('openmp', meshsize, self.mpicmd)

    def _pyop2_call(self, backend, meshsize):
        if self.profile:
            prof = '-m cProfile -o %s ' % self._path('pyop2_adv_diff.%s.cprofile' % backend)
        else:
            prof = ''
        cmd = 'python %spyop2_adv_diff.py -m %s -b %s' % (prof, self.mesh % meshsize, backend)
        return self.logged_call_with_time(cmd)

    def pyop2_seq(self, meshsize):
        """PyOP2 advection-diffusion benchmark (sequential backend)"""
        return self._pyop2_call('sequential', meshsize)

    def pyop2_openmp(self, meshsize):
        """PyOP2 advection-diffusion benchmark (OpenMP backend)"""
        return self._pyop2_call('openmp', meshsize)

    def pyop2_cuda(self, meshsize):
        """PyOP2 advection-diffusion benchmark (CUDA backend)"""
        return self._pyop2_call('cuda', meshsize)

    def __init__(self, np=1, message='', version=None, reference=None,
            extrude=False, meshsize=None, parameters=None, mpicmd=None,
            flcmd=None, profile=None):
        parameters = parameters or ['version', 'meshsize']
        # Execute for all combinations of these parameters
        self.meshsize = meshsize or ['0.000008', '0.000004', '0.000002']
        self.version = version or ['fluidity', 'fluidity_pyop2_seq', 'pyop2_seq']
        # Compare timings against this parameter value
        self.reference = ('version', reference or self.version[0])
        super(AdvDiffBenchmark, self).__init__(parameters)
        self.primed = []

        self.log("Running versions: %s" % self.version)
        self.log("Reference %s: %s" % self.reference)

        self.np=np
        self.message=message
        self.mpicmd = mpicmd or 'mpiexec' if np>1 else ''
        self.flcmd = flcmd or '${FLUIDITY_DIR}/bin/fluidity'
        self.profile = '.profile' if profile else ''
        self.meshdir = os.path.join('meshes', str(self.np)) if self.np > 1 else 'meshes'
        self.flmldir = os.path.join('flmls', str(self.np)) if self.np > 1 else 'flmls'
        if extrude:
            self.mesh = os.path.join(self.meshdir,'mesh.%s')
            self.generate_mesh = generate_meshfile
        else:
            self.mesh = os.path.join(self.meshdir,'square.%s')
            self.generate_mesh = generate_trianglefile

        self.plotstyle = dict(zip(self.version, ['k-o', 'g-s', 'r-d', 'b-^']))
        self.plotlabels = {
                'fluidity': 'Fluidity (cores: 1)',
                'fluidity_mpi': 'Fluidity (cores: %d)' % self.np,
                'fluidity_pyop2_seq': 'Fluidity-PyOP2 (backend: sequential)',
                'fluidity_pyop2_cuda': 'Fluidity-PyOP2 (backend: cuda)',
                'fluidity_pyop2_openmp': 'Fluidity-PyOP2 (backend: openmp)',
                'fluidity_pyop2_mpi': 'Fluidity-PyOP2 (backend: mpi)',
                'fluidity_pyop2_mpi_openmp': 'Fluidity-PyOP2 (backend: mpi+openmp)',
                'pyop2_seq': 'PyOP2 (backend: sequential)',
                'pyop2_openmp': 'PyOP2 (backend: OpenMP)',
                'pyop2_cuda': 'PyOP2 (backend: CUDA)',
                'dolfin': 'DOLFIN (cores: %d)' % self.np
                }

    def create_input(self, reorder):
        if not os.path.exists(self.meshdir):
            os.makedirs(self.meshdir)
        if not os.path.exists(self.flmldir):
            os.makedirs(self.flmldir)
        for s in self.meshsize:
            mesh = self.mesh % s
            # Generate mesh
            self.log("Generating mesh %s" % mesh)
            self.log(self.generate_mesh(mesh, s, capture=True, reorder=reorder, move=self.np==1))
            # Decompose the mesh if running in parallel
            if self.np > 1:
                if reorder:
                    cmd = '${HILBERT_DIR}/bin/flredecomp -i 1 -o %d reorder_mesh_0_checkpoint decomp'
                    self.logged_call([self.mpicmd, cmd % self.np])
                    for m in glob('decomp_CoordinateMesh_*'):
                        self.logged_call('mv %s %s%s' % (m, mesh, m[21:]))
                    self.logged_call("mv reorder_mesh_CoordinateMesh_0_checkpoint.ele %s.ele" % mesh)
                    self.logged_call("mv reorder_mesh_CoordinateMesh_0_checkpoint.edge %s.edge" % mesh)
                    self.logged_call("mv reorder_mesh_CoordinateMesh_0_checkpoint.node %s.node" % mesh)
                else:
                    self.logged_call('${FLUIDITY_DIR}/bin/fldecomp -m triangle -n %s %s' \
                            % (self.np, mesh))
            self.logged_call('${FENICS_DIR}/bin/dolfin-convert %s.ele %s.xml' % (mesh, os.path.join('meshes', s)))
            def write_flml(f1, f2):
                f2.write(f1.read() % {
                    'mesh': mesh,
                    'diffusivity': parameters.diffusivity,
                    'current_time': parameters.current_time,
                    'dt': parameters.dt,
                    'endtime': (parameters.endtime - parameters.dt),
                    'backend': backend
                    })
            # Generate flml
            for backend in ['sequential', 'openmp', 'cuda']:
                with open('ufl_advection_diffusion%s.flml.template' % self.profile) as f1, \
                        open(os.path.join(self.flmldir,'ufl_advection_diffusion%s.%s.%s.flml'
                             % (self.profile, backend, s)), 'w') as f2:
                    write_flml(f1, f2)
            with open('advection_diffusion.flml.template') as f1, \
                    open(os.path.join(self.flmldir,'advection_diffusion.%s.flml'%s), 'w') as f2:
                write_flml(f1, f2)

    def dry_run(self, version, meshsize):
        if version not in self.primed:
            self.log("Performing dry run of %s using mesh size %sx%s" % (version, meshsize, meshsize))
            self.__getattribute__(version)(meshsize)
            self.primed.append(version)

    def run(self, version, meshsize):
        self.dry_run(version, meshsize)
        self.log("Running %s with mesh size %s" % (version, meshsize))
        t = self.__getattribute__(version)(meshsize)
        with open(self.mesh % meshsize + '.ele') as m:
            elems = int(m.readline().split()[0])
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
            self._plot('speedup_'+fig, pl, lambda x: x+'_speedup', 'lower right',
                       'Relative speedup over Fluidity baseline', title, format)

def get_parser(desc=AdvDiffBenchmark.__doc__):
    import argparse
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-r', '--run', action='store_true', help='run benchmarks')
    parser.add_argument('-p', '--plot', action='store_true', help='plot results')
    parser.add_argument('-d', '--dump', action='store_true', help='Pickle dump')
    parser.add_argument('-c', '--csv', action='store_true', help='Dump to CSV')
    parser.add_argument('-v', '--versions', nargs='+', help='Versions to run')
    parser.add_argument('-a', '--reference', help='Version to use as reference')
    parser.add_argument('-l', '--load', help='Pickle load from file')
    parser.add_argument('-q', '--quiet', help='Only print errors and warnings')
    parser.add_argument('-e', '--extrude', action='store_true',
            help='Use extruded structured unit square mesh')
    parser.add_argument('-s', '--create-input', action='store_true',
            help='Do not generate input files')
    parser.add_argument('--reorder', action='store_true', default=True,
            help='Reorder mesh (only when called with --create-input)')
    parser.add_argument('-n', '-np', type=int, default=1,
            help='Number of MPI process (Fluidity, DOLFIN)')
    parser.add_argument('-m', '--message', default='',
            help='Message, added to the log output')
    parser.add_argument('--mpi-cmd', help='Command to run MPI processes')
    parser.add_argument('--fluidity-cmd', help='Fluidity binary to run')
    parser.add_argument('--python-profile', action='store_true',
            help='Create a cProfile of the Python runs')
    return parser

def main(args):
    logging.getLogger().setLevel(logging.WARN if args.quiet else logging.INFO)
    b = AdvDiffBenchmark(args.n, args.message, args.versions, args.reference,
            args.extrude, mpicmd=args.mpi_cmd, flcmd=args.fluidity_cmd, profile=args.python_profile)
    if args.load and os.path.exists(args.load):
        b.load(args.load)
    if args.create_input:
        b.create_input(args.reorder)
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

if __name__ == '__main__':
    main(get_parser().parse_args())
