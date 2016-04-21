import os
import time

from pyBadlands.remote import RemoteModel
from pyBadlands.model import Model

start_time = time.time()
model = Model()
model.load_xml('bench/input.xml')
model.run_to_time(4000)
print 'singlethread finished in %s seconds' % (time.time() - start_time)

for ncpus in range(1, 5):
    start_time = time.time()
    remote = RemoteModel(maxcpus=ncpus)
    # match the child's cwd to ours
    remote._view.execute('import os; os.chdir("%s")' % os.getcwd())
    remote.load_xml('bench/input.xml')
    remote.run_to_time(4000)
    print 'mpi (%d cores) finished in %s seconds' % (ncpus, time.time() - start_time)
