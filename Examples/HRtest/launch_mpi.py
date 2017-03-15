#!/usr/bin/env python

# This script shows how to run a model on an MPI cluster. It runs the HRtest
# example.
#
# To run this on an MPI cluster, use something like:
#
#     mpiexec -n 4 python launch_mpi.py 1000

import os
import re
import sys
import time

from pyBadlands.model import Model

base_path = '.'
xml_name = 'hrtest.xml'

run_years = int(sys.argv[1])

# Modify base_path to run from another directory 
# Change into HRtest directory if required...
os.chdir(base_path)

start_time = time.time()
model = Model()

print 'loading %s' % xml_name

print model.load_xml(xml_name)

print model.run_to_time(run_years)
print 'run to %s years finished in %s seconds' % (run_years, time.time() - start_time)