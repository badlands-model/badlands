import os
import time

from pyBadlands.remote import RemoteModel
from pyBadlands.model import Model

start_time = time.time()
model = Model()
model.load_xml('bench/input.xml')
model.run_to_time(4000)
print 'singlethread finished in %s seconds' % (time.time() - start_time)

start_time = time.time()
remote = RemoteModel()
# match the child's cwd to ours
remote._view.execute('import os; os.chdir("%s")' % os.getcwd())
remote.load_xml('bench/input.xml')
print 'load time is %s' % time.time()
remote.run_to_time(4000)
print 'run time is %s' % time.time()
print 'mpi finished in %s seconds' % (time.time() - start_time)
