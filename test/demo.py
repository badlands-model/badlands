from pyBadlands.model import Model as badlandsModel

# from ipyparallel import Client

# Connect to the MPI cluster and configure our view
# client = Client(profile='mpi')

# initialise model
model = badlandsModel()
model.load_xml('data/demo-input.xml')

model.run_to_time(100)  # run model for 100 years
