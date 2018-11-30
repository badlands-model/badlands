#!/bin/bash

ipcluster start --profile=mpi &
jupyter notebook --ip=0.0.0.0 --no-browser --NotebookApp.token='' --allow-root
