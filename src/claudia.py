
# Imports

# #+srcname: claudia-imports

import numpy
import argparse
import string

# ** The SmartDict class

# + Taken from the excellent [[http://code.activestate.com/recipes/577590-dictionary-whos-keys-act-like-attributes-as-well/][Python Recipe]] by [[http://code.activestate.com/recipes/users/4174115/][Sunjay Varma]]
# + I will use this as a base class 

# #+srcname: claudia-smartdict

class SmartDict(dict):
    """
    Combines the best features of a class and a dict
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(e)
    def __setattr__(self, name, value):
        self[name] = value

# ** The class for a Cloudy model

# #+srcname: claudia-model-class

class CloudyModel(SmartDict):
    """
    A single Cloudy model
    """
    indir, outdir = "in", "out"
    insuff, outsuff = ".in", ".out"
    def __init__(modelname):
        self.infilepath = os.join(indir, modelname + insuff)
        with open(self.infilepath) as f:
            self._inscript = f.read() 
        savecommands = find_save_commands(self._inscript)

# *** List of possibilities for cloudy save files

# + Taken from Hazy1 C10 version 2011/08/14
# + This is nowhere near exhaustive
# + These are checked in turn, so more specific types should come first. 

# #+srcname: claudia-types-of-cloudy-save-files

SAVETYPES = [
    "cooling",
    "element hydrogen",
    "element helium",
    "element carbon",
    "element nitrogen",
    "element oxygen",
    "element sulfur",
    "element silicon",
    "element iron",
    "heating",
    "lines emissivity",
    "overview",
    "PDR",
    "physical conditions",
    "pressure",
    "source function, spectrum",
    "source function, depth",
    ]

# *** Find basic info about the run
#     :LOGBOOK:
#     CLOCK: [2011-08-20 Sat 18:24]--[2011-08-21 Sun 00:04] =>  5:40
#     :END:

# #+srcname: claudia-input-parse-basic-info



# *** Find which save files were written
#     :LOGBOOK:
#     - Note taken on [2011-08-20 Sat 18:21] \\
#       OK, this is just about working now, time to move on
#     - Note taken on [2011-08-20 Sat 14:16] \\
#       Not sure what we were doing here? What was the use-case of the cut_out function.
#     CLOCK: [2011-08-20 Sat 14:16]--[2011-08-20 Sat 18:24] =>  4:08
#     CLOCK: [2011-06-28 Tue 13:14]--[2011-06-28 Tue 13:16] =>  0:02
#     CLOCK: [2011-06-27 Mon 23:46]--[2011-06-27 Mon 23:46] =>  0:00
#     :END:

# This originally seemed like a job for regular expressions, but that quickly got out of hand. 

# Instead of allowing any type of save file, we use a finite list =SAVETYPES= since that makes the parsing much simpler. The only problem is that Cloudy allows the names to be abbreviated to four letters. 

# #+srcname: claudia-get-list-of-save-files

def find_save_commands(s):
    """
    Find all save commands in a Cloudy input file and return a list of [type, file] pairs

    >>> find_save_commands('save heating last ".heat"\\nsave cooling last ".cool"')
    [('heating', '.heat'), ('cooling', '.cool')]
    """
    save_commands = [] 
    for line in s.split("\n"):
        found = find_single_save_command(line)
        if found: save_commands.append(found)
    return save_commands or None
    

def find_single_save_command(line):
    """
    Parse single line of a Cloudy input file, looking for a save command

    It should work both with C08-style (punch) and C10-style (save) commands:

    >>> find_single_save_command('save overview last ".ovr"')
    ('overview', '.ovr')
    >>> find_single_save_command('PUNCH LAST OVERVIEW ".ovr"')
    ('overview', '.ovr')
    >>> find_single_save_command('save over no buffering, last, file=".ovr"')
    ('overview', '.ovr')
    >>> find_single_save_command('save madeupname file=".xyz"')
    (None, '.xyz')
    >>> find_single_save_command('this is not the right command')

    Note that the last command prints nothing since it returns None
   
    """
    line = line.lower()
    if line.startswith("save") or line.startswith("punch"):
        assert '"' in line or "'" in line, "No filename given in save/punch command"
        line = cut_out(line, "save")
        line = cut_out(line, "punch")
        if "last" in line:
            line = cut_out(line, "last")
        if '"' in line:
            delim = '"'
        elif "'" in line:
            delim = "'"
        firstpart, savefile = line.split(delim)[:2]
        for savetype in SAVETYPES:
            if look4stringinline(savetype, firstpart):
                return savetype, savefile
        # failed to find anything
        return None, savefile
    else:
        return None

# *** Utility functions for input parsing 
# #+srcname: claudia-input-parse-utilities

def cut_out(s, phrase):
    """
    Returns the input string <s> but with all occurrences of <phrase> deleted

    <phrase> should be one or more words, separated by whitespace. Effort is made
    to preserve one space between words, which makes it better than s.replace(phrase, '')

    >>> s = 'the quick brown fox, which is the brownest ever, jumped over the lazy dog'
    >>> cut_out(s, 'the')
    'quick brown fox, which is brownest ever, jumped over lazy dog'
    >>> s.replace('the', '')
    ' quick brown fox, which is  brownest ever, jumped over  lazy dog'

    Note the extra spaces in the s.replace version
    """
    return ' '.join(map(string.strip, s.split(phrase))).strip()

def look4stringinline(string, line):
    """
    Look for string in line, only comparing the first 4 characters of each word

    This is because cloudy does the same.

    Case should not matter: 
    >>> look4stringinline('punch pressure', 'PUNC FINAL PRES')
    True

    And it is OK to have strings with less than 4 characters:
    >>> look4stringinline('PDR', 'save pdr')
    True

    And here is an example that should fail:
    >>> look4stringinline('save whatever', 'save foobar')
    False

    """
    words = string.split()
    for word in words:
        if len(word) > 4: word = word[:4] 
        if not word.upper() in line.upper():
            return False
    return True
