from pyBadlands.model import Model as badlandsModel

# from ipyparallel import Client

# Connect to the MPI cluster and configure our view
# client = Client(profile='mpi')

outDir = 'output/'

# initialise model
model = badlandsModel()

# use default parameters (see model.__init__)

output_dir = 'output'  # FIXME: why do we need this at load time?
model.load_dem('data/regularMR.csv')  # load input file, generate TIN

# configure rainfall pattern
rainVal = 1.  # Precipitation rate [m/a]
model.rain.fill(rainVal)
model.rain[:model.recGrid.boundsPt] = 0.

# run 100kyears in 1k year intervals
# for t in xrange(0, 100000, 1000):
for t in xrange(1, 31, 5):
    print('run to t=%d' % t)
    model.compute_flow(tEnd=t)

    step = t / 5  # TODO: confirm what this means - increasing integer for vis?
    model.write_output(outDir=output_dir, step=step)  # write HDF5 output
