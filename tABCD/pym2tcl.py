#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '..'))
import msutils

msutils.make_dir("output/tcl")

for i in range(10000):
    msutils.pym2tcl("output/models/configuration.%d.pym" % i,
                    "output/tcl/configuration.%d.tcl" % i)
