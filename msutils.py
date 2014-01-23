"""Utilities used by most scripts in the benchmark"""

import os
import sys
import re

def make_dir(dirname):
    """Make a directory (and parents) if it doesn't exist."""
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def pym2tcl(pym_name, tcl_name):
    """Generate a tcl file from a pym file."""
    if not os.path.exists(pym_name):
        return
    InputFileLines = open(pym_name).read().splitlines()

    f1 = open(tcl_name, 'w')
    f1.write('draw color blue' + '\n')
    for line in InputFileLines:
        if line.find('SPHERE') >= 0:
            seq = line[:]

            Newseq = seq.replace(",", " ")
            l = re.findall("(\S+)", Newseq)
            seqX = str(float(l[1]))
            seqY = str(float(l[2]))
            seqZ = str(float(l[3]))
            seqRad = str(float(l[4]))
            f1.write('draw sphere { ' + seqX + ' ' + seqY + ' ' +
                     seqZ + ' ' + '} radius ' + seqRad +
                     ' resolution 300' + '\n')
