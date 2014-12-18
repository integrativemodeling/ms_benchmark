#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '..', '..'))
import msutils

msutils.make_dir("output/mfj")

for i in range(10000):
    msutils.tcl2mfj("output/tcl/configuration.%d.tcl" % i,
                    "output/mfj/configuration.%d.mfj" % i, 5)
