from ipyparallel import Client


def relog():
    ''' For debugging, redirect the individual node stdout/stderr to a file '''
    import os
    pid = os.getpid()
    logfile = open('/tmp/model-%s.txt' % pid, 'w')
    logfile.write('--- I am PID %s\n' % pid)

    import sys
    sys.stdout = logfile
    sys.stderr = logfile


class RemoteModel(object):
    """
    Wrapper to allow Model to run on an MPI cluster while hiding most of the
    MPI details.

    The public interface is identical to Model, so you should be able to
    substitute one for the other as desired..

    We use MPI in a master-slave architecture. This object runs on the master
    and handles all communication with the slaves. The slaves run the native
    Model object. In this way, the calling code does not need to be aware of the
    underlying parallelisation details.
    """

    def __init__(self, profile='mpi'):
        client = Client(profile=profile)
        self._view = client[:]
        self._view.block = True

        self._view.execute('import os')
        self._view.execute('cwd = os.getcwd()')

        self._view.execute('from pyBadlands.model import Model')
        self._view.execute('model = Model()')

        # Uncomment this to enable node debug logging to /tmp
        # self._view.apply(relog)

    def load_xml(self, filename, verbose=False):
        try:
            self._view.execute('model.load_xml(filename="%s", verbose=%s)' % (filename, verbose))
        except Exception, e:
            import pdb; pdb.set_trace()

    def run_to_time(self, tEnd):
        try:
            self._view.execute('model.run_to_time(%s)' % tEnd)
        except Exception, e:
            import pdb; pdb.set_trace()
