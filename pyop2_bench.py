import logging
import os
import csv
import cPickle
from collections import defaultdict
from datetime import datetime
import matplotlib as mpl
mpl.use("Agg")
import numpy as np
import pylab
import subprocess

from benchrun import Benchmark


class PyOP2Benchmark(Benchmark):
    """Abstract base class for PyOP2 benchmarks."""

    def __init__(self, parameters):
        super(PyOP2Benchmark, self).__init__(parameters)
        self.plotdata = defaultdict(list)
        self.timestamp = datetime.today().isoformat().replace(':', '')
        logging.basicConfig(format='%(message)s')

    def time_all(self):
        self.logged_call('uname -a')
        self.log('Timestamp: %s' % self.timestamp)
        self.log('=== Start Benchmark run at %s ===\n' % datetime.today())
        super(PyOP2Benchmark, self).time_all()
        self.log('=== Finish Benchmark run at %s ===\n' % datetime.today())

    def _path(self, filename):
        if not os.path.exists(self.timestamp):
            os.makedirs(self.timestamp)
        return os.path.join(self.timestamp, filename)

    def write_csv(self):
        with open(self._path('results.csv'), 'wb') as f:
            w = csv.writer(f)
            w.writerows(self.results)

    def log(self, msg=''):
        with open(self._path('benchmark.log'), 'a') as f:
            f.write(msg + '\n')
        logging.info(msg)
        return msg

    def logged_call(self, call, cwd=None):
        if not isinstance(call, str):
            call = ' '.join(call)
        return self.log(subprocess.check_output(call, stderr=subprocess.STDOUT, shell=True, cwd=cwd))

    def logged_call_with_time(self, call, env=None):
        if env is None:
            env = os.environ
        msg = self.logged_call(call)
        for line in msg.split('\n'):
            if line.find('/fluidity ::') != -1:
                    time = float(line.split(' ')[2])
        return time

    def flufl_call_with_time(self, call):
        msg = self.logged_call(call)
        for line in msg.split('\n'):
            if line.find('UFL ::') != -1:
                    time = float(line.split(' ')[2])
        return time

    def dump(self):
        with open(self._path('results.pickle'), 'wb') as f:
            cPickle.dump(self.__dict__, f, cPickle.HIGHEST_PROTOCOL)

    def load(self, filename):
        with open(filename, 'rb') as f:
            self.__dict__.update(cPickle.load(f))

    def compute_speedup(self):
        for v in self.version:
            for this, base in zip(self.plotdata[v], self.plotdata[self.reference[1]]):
                self.plotdata[v + '_speedup'].append(base / this)

    def _plot(self, figname, plot, col, legend_pos, ylabel, title,
              format='svg', write_script=True):
        fig = pylab.figure(figname, figsize=(8, 6), dpi=300)
        for v in self.version:
            style, lw, label = self.plotmeta[v]
            plot(self.plotdata['elements'], self.plotdata[col(v)],
                 style, lw=lw, label=label)
        pylab.legend(loc=legend_pos)
        pylab.xlabel('Number of elements in the mesh')
        pylab.ylabel(ylabel)
        pylab.title(title)
        pylab.grid()
        if not format:
            pylab.show()
        else:
            for fmt in format.split(','):
                pylab.savefig(self._path('%s.%s' % (figname, fmt)),
                              orientation='landscape', format=fmt, transparent=True)
        pylab.close(fig)

        if write_script:
            np.save(self._path('elements.npy'), self.plotdata['elements'])
            code = "pylab." + plot.func_name + \
                   "(np.load('elements.npy'), np.load('%s.npy'), '%s', lw=%d, label='%s')\n"
            plots = ''
            for v in self.version:
                np.save(self._path(col(v) + '.npy'), self.plotdata[col(v)])
                style, lw, label = self.plotmeta[v]
                plots += code % (col(v), style, lw, label)
            if not format:
                savefig = "pylab.show()"
                mplimport = ""
            else:
                code = "pylab.savefig('%s.%s', orientation='landscape', format='%s', transparent=True)"
                savefig = '\n'.join([code % (figname, fmt, fmt) for fmt in format.split(',')])
                mplimport = "import matplotlib\nmatplotlib.use('Agg')"
            with open(self._path(figname + ".plot.py"), 'w') as f:
                f.write("""import numpy as np
%(mplimport)s
import pylab

fig = pylab.figure('%(figname)s', figsize=(8, 6), dpi=300)
%(plots)s
pylab.legend(loc='%(legend_pos)s')
pylab.xlabel('Number of elements in the mesh')
pylab.ylabel('%(ylabel)s')
pylab.title('%(title)s')
pylab.grid()
%(savefig)s
""" % locals())
