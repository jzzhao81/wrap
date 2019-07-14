#!/usr/bin/env python
# coding=UTF-8
'''
@Author: Jianzhou Zhao
@Date: 2019-07-14 10:27:09
@LastEditors: Jianzhou Zhao
@LastEditTime: 2019-07-14 11:10:33
@Description: This code prepares the input of DMFT from a Wannier tight-binding (TB).
'''

import os
from collections import defaultdict

import mpiutils as mpi
import numpy as np


class Model():

    def __init__(self, nwan=None, nrpt=None, deg=None, ham=None):
        self.nwan = nwan
        self.nrpt = nrpt
        self.deg = deg
        self.ham = ham
        return

    @classmethod
    def from_hr(cls, hr_file):

        if mpi.rank == 0:
            with open(hr_file, 'r') as f:
                next(f)  # skip title line
                nwan = int(next(f))
                nrpt = int(next(f))

                # read degeneration of R points
                deg = []
                for _, line in zip(range(int(np.ceil(nrpt / 15))), f):
                    deg.extend(int(i) for i in line.split())

                assert len(deg) == nrpt

                # read rpt and hamr
                raw_list = [line.split() for line in f]

                nwan_square = nwan**2
                hr_list = defaultdict(list)
                for num, line in enumerate(raw_list):
                    Rvec = tuple(map(int, line[:3]))
                    hr_list[Rvec].append(
                        (float(line[5]) + 1j * float(line[6])) / deg[num // nwan_square])

            for key, array in hr_list.items():
                hr_list[key] = [array[iwan::nwan] for iwan in range(nwan)]
        else:
            nwan = 0
            nrpt = 0
            deg = None
            hr_list = None

        nwan = mpi.bcast(nwan, root=0)
        nrpt = mpi.bcast(nrpt, root=0)
        deg = mpi.bcast(deg, root=0)
        hr_list = mpi.bcast(hr_list, root=0)

        return cls(nwan=nwan, nrpt=nrpt, deg=deg, ham=hr_list)

    def get_bulk_Hk(self, Kvec):
        return np.sum([np.array(Hvec, dtype=np.complex) * np.exp(2j * np.pi * np.dot(Rvec, Kvec))
                       for Rvec, Hvec in self.ham.items()], axis=0)

    def get_eigval(self, Kvec):
        import scipy.linalg as la
        return la.eigvalsh(self.get_bulk_Hk(Kvec))

    def get_eigvals(self, klist):
        kmpi = mpi.devide_array(klist, root=0)
        empi = np.array([self.get_eigval(kvec) for kvec in kmpi])
        eigs = np.zeros((len(klist), self.nwan), dtype=np.float)
        mpi.Gatherv(empi, (eigs, len(empi.flatten())), root=0)
        return eigs.T


class Kpoint():

    def __init__(self, nktot=None, kfrac=None, kreal=None):
        self.nktot = nktot
        self.kfrac = kfrac
        self.kreal = kreal
        return

    @staticmethod
    def k_length(K_list):
        """
        Return k point distance from previous one
        Exclude the first point
        """
        return np.array([np.linalg.norm(K1 - K0) for K0, K1 in zip(K_list, K_list[1:])])

    @classmethod
    def k_length_accumulate(cls, K_list):
        """
        Return k point distance from first one
        Include the first pointself.

        This code should be written beautifully!

        """
        dist = [0.0]
        dist.extend(cls.k_length(K_list))
        return [sum(dist[:ii])+dist[ii] for ii in range(len(dist))]

    @classmethod
    def from_high_symmetry_points(cls, khsym_frac, bzvec=None, kwant=None):
        """
        Return a K list along high symmetry lines given.
        """
        if bzvec is None:
            bzvec = np.eye(3)

        khsym_real = np.array([np.dot(kvec, bzvec) for kvec in khsym_frac])
        klen = np.sum(cls.k_length(khsym_real))

        if kwant is None or kwant < len(khsym_frac):
            kwant = len(khsym_frac) - 1
        kstep = klen / kwant

        kline_real = []
        for k0, k1 in zip(khsym_real, khsym_real[1:]):
            nk = int(np.ceil(np.linalg.norm(k1 - k0) / kstep))
            kfrac = (np.array(k1) - np.array(k0)) / nk
            kline_real.extend(list(k0 + ik * kfrac) for ik in range(nk))
        kline_real.append(list(khsym_real[-1]))

        kline_frac = []
        for kpt in kline_real:
            kline_frac.append(list(np.dot(kpt, np.linalg.inv(bzvec))))

        return cls(nktot=len(kline_frac), kfrac=np.array(kline_frac), kreal=np.array(kline_real))

    @classmethod
    def gen_uniform_kpoints(cls, nk, bzvec=None):
        """
        Return a uniform K point grid.
        """
        if bzvec is None:
            bzvec = np.eye(3)
        nktot = np.prod(nk)
        kfrac = [[ikx/nk[0], iky/nk[1], ikz/nk[2]] for ikz in range(nk[2]) for iky in range(nk[1]) for ikx in range(nk[0])]
        kreal = [np.dot(kpt, bzvec) for kpt in kfrac]
        return cls(nktot=nktot, kfrac=np.array(kfrac), kreal=np.array(kreal))


def gen_proj_from_string(w90file='wannier90_hr.dat', nkpt=(4, 4, 4), pjorb=None, occ=0.0):

    pass


if __name__ == "__main__":

    gen_proj_from_string()
