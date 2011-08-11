
import numpy
import argparse
import string

def find_save_commands(s):
    """
    Find all save commands in a Cloudy input file and return a list of [type, file] pairs

    >>> find_save_commands('save overview last ".ovr"\\nsave pressures last ".pre"')
    [['overview', '.ovr'], ['pressures', '.pre']]
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
    ['overview', '.ovr']
    >>> find_single_save_command('punch last overview ".ovr"')
    ['overview', '.ovr']
    """
    if line.startswith("save") or line.startswith("punch"):
        pass
    else:
        return None

def cut_out(s, phrase):
    """
    Returns the input string <s> but with all occurences of <phrase> deleted

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
