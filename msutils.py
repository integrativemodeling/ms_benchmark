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

def tcl2mfj(tcl_name, mfj_name, num_spheres):
    """Generate an mfj file from a tcl file."""
    if not os.path.exists(tcl_name):
        return
    Nsphs = str(num_spheres)
    lines = open(tcl_name).read().splitlines()

    f1 = open(mfj_name, 'w')
    f1.write(os.path.basename(mfj_name) + '\n')
    f1.write('1' + '\n')
    f1.write(Nsphs + '\n')
    f1.write('ang' + '\n')
    f1.write('none' + '\n')
    f1.write('1.0000' + '\n')
    for line in lines:
        if line[0:11] == 'draw sphere':
            l = re.findall("(\S+)", line)
            x = float(l[3])
            y = float(l[4])
            z = float(l[5])
            r = float(l[8])
            f1.write("%s%8.3f%s%8.3f%s%8.3f%7.1f%s" %
                     ('    ', x, '      ', y, '      ',
                      z, r, '      .00000') + '\n')
