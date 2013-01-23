import logging
import os
import csv
import pickle
from collections import defaultdict
from datetime import datetime
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
            f.write(msg+'\n')
        logging.info(msg)

    def logged_call(self, call):
        self.log(subprocess.check_output(call, stderr=subprocess.STDOUT, shell=True))

    def dump(self):
        with open(self._path('results.pickle'), 'wb') as f:
            pickle.dump(self.__dict__, f)

    def load(self, filename):
        with open(filename, 'rb') as f:
            self.__dict__.update(pickle.load(f))

    def compute_speedup(self):
        for v in self.version:
            for this, base in zip(self.plotdata[v], self.plotdata[self.reference[1]]):
                self.plotdata[v+'_speedup'].append(base/this)

    def _plot(self, fig, plot, col, legend_pos, ylabel, title, format='svg'):
        f = pylab.figure(fig, figsize=(8, 6), dpi=300)
        for v in self.version:
            plot(self.plotdata['elements'], self.plotdata[col(v)], self.plotstyle[v], lw=2, label=self.plotlabels[v])
        pylab.legend(loc=legend_pos)
        pylab.xlabel('Number of elements in the mesh')
        pylab.ylabel(ylabel)
        pylab.title(title)
        pylab.grid()
        if not format:
            pylab.show()
        else:
            for fmt in format.split(','):
                pylab.savefig(self._path('%s.%s' % (fig,fmt)), orientation='landscape', format=fmt, transparent=True)
        pylab.close(f)
