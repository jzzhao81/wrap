#!/usr/bin/env python
# coding=UTF-8
'''
@Author: Jianzhou Zhao
@LastEditors: Jianzhou Zhao
@Description: This is a wrapped code for the mpi4py
@Date: 2019-04-19 22:00:22
@LastEditTime: 2019-07-12 15:02:16
'''

import sys
import warnings

rank = 0
size = 1
comm = None
world = None
rank0 = True

# Try to setup MPI and get the comm, rank and size.
# If not they should end up as rank=0, size=1.
try:

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    double = MPI.DOUBLE

    if comm is not None and size > 1:
        print("Starting MPI rank=%i [size=%i]" % (rank, size))

    rank0 = True if rank == 0 else False

except ImportError:

    warnings.warn("Warning: mpi4py not installed.")


def bcast(obj, root=0, comm=comm):
    if comm is not None and comm.size > 1:
        return comm.bcast(obj, root=root)
    else:
        return obj


def reduce(obj, root=0, op=None, comm=comm):
    if comm is not None and comm.size > 1:
        return comm.reduce(obj, root=root, op=(op or MPI.SUM))
    else:
        return obj


def allreduce(obj, op=None, comm=comm):
    if comm is not None and comm.size > 1:
        return comm.allreduce(obj, op=(op or MPI.SUM))
    else:
        return obj


def gather(obj, root, comm=comm):
    if comm is not None and comm.size > 1:
        return comm.gather(obj, root=root)
    else:
        return obj


def Reduce(sendobj, recvobj, root=0, op=None, comm=comm):
    if comm is not None and comm.size > 1:
        comm.Reduce(sendobj, recvobj, root=root, op=(op or MPI.SUM))
    else:
        recvobj = sendobj


def Allreduce(sendobj, recvobj, op=None, comm=comm):
    if comm is not None and comm.size > 1:
        comm.Allreduce(sendobj, recvobj, op=(op or MPI.SUM))
    else:
        recvobj = sendobj


def Gather(sendobj, recvobj, root, comm=comm):
    if comm is not None and comm.size > 1:
        comm.Gather(sendobj, recvobj, root=root)
    else:
        recvobj = sendobj


def Gatherv(sendobj, recvobj, root, comm=comm):
    comm.Gatherv(sendbuf=sendobj, recvbuf=recvobj, root=root)


def Scatterv(sendobj, recvobj, root, comm=comm):
    if comm is not None and comm.size > 1:
        comm.Scatterv(sendobj, recvobj, root=root)
    else:
        recvobj = sendobj


def devide_array(sendobj, root):
    """
    Devide sendobj into size part.
    """
    from numpy import zeros
    from functools import reduce

    total_count = len(sendobj)
    " For array more than 1 dimension, we only devide the first axis "
    coeff = reduce(lambda x, y: x*y, sendobj.shape[1:])
    count = [coeff*(total_count//size+1) if irank <= total_count % size-1
             else coeff*(total_count//size) for irank in range(size)]
    disp = [sum(count[:irank]) for irank in range(size)]
    recvobj = zeros(count[rank])
    comm.Scatterv([sendobj, count, disp, MPI.DOUBLE], recvobj, root=root)
    return recvobj.reshape(-1, *sendobj.shape[1:])
