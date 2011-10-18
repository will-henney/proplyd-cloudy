"""
Define a lightweight Model() class for running Cloudy models (optionally in parallel)

"""


import os
import subprocess
import multiprocessing
import sys

class _InputScript(str):
    """
    Input script for a single Cloudy model
    """
    def __init__(self):
        pass
    
    
    
class Model(object):
    """
    A single cloudy model.

    Public methods are:
    write() : used to add lines to the input script
    run() : used to run the model (after first writing input script to a file)
    
    """
    indir = "in"
    outdir = "out"
    cloudy_cmd = ["time", "cloudy.exe"] 
    verbose = True

    def __init__(self, prefix, **kwargs):
        self.prefix = prefix
        # Any optional keywords get set as attributes
        self.__dict__.update(kwargs)
        self._create_dirs_as_necessary()
        self.infilepath = '%s/%s.in' % (self.indir, self.prefix)
        self.outfilepath = '%s/%s.out' % (self.outdir, self.prefix)
        self.input = _InputScript()
        # The following is necessary so that any changes to the class variable are remembered
        # when the instances are pickled for multiprocessing
        self.cloudy_cmd = Model.cloudy_cmd
        self._writeheader()

    def _create_dirs_as_necessary(self):
        
        

    def run(self):
        self._write_input_script_to_file()
        status = subprocess.Popen(self.cloudy_cmd, 
                                  stdin=file(self.infilepath, "r"),
                                  stdout=file(self.outfilepath, "w"),
                                  stderr=subprocess.STDOUT,
                                  cwd=self.outdir
                                  ).wait()
        if self.verbose:
            try:
                print multiprocessing.current_process().name, " finished with ", self.prefix
            except:
                pass
        return status

    def __call__(self):
        """
        This is needed because bound methods such as run() cannot be
        pickled and so cannot be passed to multiprocessing
        pools. Instead, we define this __call__ method so that the
        class instance is callable, then we send the instance itself
        (which *is* pickleable) directly to the pool. 
        
        """
        return self.run()

    def _write_input_script_to_file(self):
        with file(self.infilepath, "w") as f:
            f.write(self.input)

    
    def write(self, s):
        self.input += s
        
    def _writeheader(self):
        self.write("* Cloudy input file written by %s\n" % (sys.argv[0]))
        self.write("title %s\n" % (self.prefix))
        self.write('set save prefix "%s"\n' % (self.prefix))
        
