#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import glob

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def test_aabc(self):
        """Test AABC example"""
        os.chdir(os.path.join(TOPDIR, 'tAABC'))
        self.run_scripts(min_models=9900, min_threshold=2500)

    def test_abcd(self):
        """Test ABCD example"""
        os.chdir(os.path.join(TOPDIR, 'tABCD'))
        self.run_scripts(min_models=9900, min_threshold=2500, threshold=2)

    def test_abcde(self):
        """Test ABCDE example"""
        os.chdir(os.path.join(TOPDIR, 'tABCDE'))
        self.run_scripts(min_models=9900, min_threshold=2500, threshold=4)

    def test_aabbc(self):
        """Test AABBC example"""
        os.chdir(os.path.join(TOPDIR, 'tAABBC'))
        self.run_scripts(min_models=9900, min_threshold=2500, threshold=5)

    def test_11mer(self):
        """Test 11mer example"""
        os.chdir(os.path.join(TOPDIR, 't-11mer'))
        self.run_scripts(min_models=4900, min_threshold=1000, threshold=15)

    def test_clamp_7mer(self):
        """Test clamp loader 7mer example"""
        os.chdir(os.path.join(TOPDIR, 'clamp-loader', 'clamp-loader-7mer'))
        self.run_scripts(min_models=9900, min_threshold=2500, threshold=11)

    def test_clamp_5mer(self):
        """Test clamp loader 5mer example"""
        os.chdir(os.path.join(TOPDIR, 'clamp-loader', 'clamp-loader-5mer'))
        self.run_scripts(min_models=9900, min_threshold=2500, threshold=6)

    def run_scripts(self, min_models, min_threshold, threshold=1):
        # Test modeling script
        p = subprocess.check_call(['./ms_cg.py'])
        # Make sure models and clusters were produced
        num = len(glob.glob('output/models/configuration.*.pym'))
        self.assertTrue(num > min_models, "Only %d models were produced" % num)
        self.assertTrue(len(glob.glob('output/models/cluster.*.pym')), 10)
        # Make mfj files
        p = subprocess.check_call(['./pym2tcl.py'])
        p = subprocess.check_call(['./tcl2mfj.py'])
        # Check the models against known connectivity
        p = subprocess.check_call(['./score_connect.py'])
        f = open('output/score-connect.txt').readlines()[1]
        val = [float(x) for x in f.rstrip('\r\n').split(' ')]
        num = len(val)
        self.assertTrue(num > min_models, "Only %d models were checked" % num)
        num_threshold = {}
        for t in range(1, 30):
            num_threshold[t] = len([x for x in val if x < t])
        self.assertTrue(num_threshold[threshold] > min_threshold,
                        "Only %d models were below %d (%s)" % \
                        (num_threshold[threshold], threshold,
                         ", ".join("%d < %d" % (num_threshold[x], x) \
                                                for x in range(1, 30))))

if __name__ == '__main__':
    unittest.main()
