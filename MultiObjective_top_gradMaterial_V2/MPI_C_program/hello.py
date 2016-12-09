#hello.py
from mpi4py import MPI

import mpi4py
mpi4py.rc.recv_mprobe = False

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print ("hello world from process ", rank)

