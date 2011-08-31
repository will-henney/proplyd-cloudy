import os
import subprocess

class InputScript(str):
    """
    Input script for a single Cloudy model
    """
    def __init__(self):
        pass
    
    
    
class Model(object):
    """
    A single cloudy model to run
    """

    indir = "in"
    outdir = "out"
    cloudy_cmd = ["time", "cloudy.exe"] 


    def __init__(self, prefix):
        self.prefix = prefix
        self.infilepath = '%s/%s.in' % (self.indir, self.prefix)
        self.outfilepath = '%s/%s.out' % (self.outdir, self.prefix)
        self.input = InputScript()

    def run(self):
        self._write_input_script_to_file()
        return subprocess.Popen(self.cloudy_cmd, 
                                stdin=file(self.infilepath, "r")
                                stdout=file(self.outfilepath, "w")
                                stderr=subprocess.STDOUT,
                                cwd=self.outdir
                                ).wait()

    def _write_input_script_to_file(self):
        with file(self.infilepath, "w") as f:
            f.write(self.input)

    
    def write(self, s):
        self.input += s
        
