#!/usr/bin/env python
# coding=UTF-8
'''
@Author: Jianzhou Zhao
@Date: 2019-07-14 10:27:09
@LastEditors: Jianzhou Zhao
@LastEditTime: 2019-07-15 17:12:51
@Description: This code prepares the input of DMFT from a Wannier tight-binding (TB).
'''

import os
import itertools
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
    def from_hr(cls, hr_file, spincopy=False):

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
                array = np.array(array).reshape(nwan, nwan)
                if spincopy:
                    hr_list[key] = np.kron(np.eye(2), array)
                else:
                    hr_list[key] = array

            nwan = nwan*2 if spincopy else nwan

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

    def get_eigvec(self, Kvec):
        import scipy.linalg as la
        return la.eigh(self.get_bulk_Hk(Kvec))

    def get_eigvals(self, klist):
        kmpi = mpi.devide_array(klist, root=0)
        empi = np.array([self.get_eigval(kvec) for kvec in kmpi])
        eigs = np.zeros((len(klist), self.nwan), dtype=np.float)
        mpi.Gatherv(empi, (eigs, len(empi.flatten())), root=0)
        return eigs.T

    def get_eigvecs(self, klist):
        kmpi = mpi.devide_array(klist, root=0)
        empi = []
        vmpi = []
        for kvec in kmpi:
            etmp, vtmp = self.get_eigvec(kvec)
            empi.append(etmp)
            vmpi.append(vtmp)
        empi = np.array(empi)
        vmpi = np.array(vmpi)
        eigs = np.zeros((len(klist), self.nwan), dtype=np.float)
        vecs = np.zeros((len(klist), self.nwan, self.nwan), dtype=np.complex)
        mpi.Gatherv(empi, (eigs, len(empi.flatten())), root=0)
        mpi.Gatherv(vmpi, (vecs, len(vmpi.flatten())), root=0)
        return eigs, vecs


class Kpoint():

    def __init__(self, nkpt=None, frac=None, real=None):
        self.nkpt = nkpt
        self.frac = frac
        self.real = real
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
            frac = (np.array(k1) - np.array(k0)) / nk
            kline_real.extend(list(k0 + ik * frac) for ik in range(nk))
        kline_real.append(list(khsym_real[-1]))

        kline_frac = []
        for kpt in kline_real:
            kline_frac.append(list(np.dot(kpt, np.linalg.inv(bzvec))))

        return cls(nkpt=len(kline_frac), frac=np.array(kline_frac), real=np.array(kline_real))

    @classmethod
    def gen_uniform_kpoints(cls, kgrid, bzvec=None):
        """
        Return a uniform K point grid.
        """
        if bzvec is None:
            bzvec = np.eye(3)
        nkpt = np.prod(kgrid)
        frac = [[ikx/kgrid[0], iky/kgrid[1], ikz/kgrid[2]] for ikz in range(kgrid[2]) for iky in range(kgrid[1]) for ikx in range(kgrid[0])]
        real = [np.dot(kpt, bzvec) for kpt in frac]
        return cls(nkpt=nkpt, frac=np.array(frac), real=np.array(real))


def gen_proj_from_string(w90file='wannier90_hr.dat', kgrid=(2, 2, 2), pjorb=None, spincopy=False):

    model = Model.from_hr(w90file, spincopy=spincopy)
    kpts = Kpoint.gen_uniform_kpoints(kgrid=kgrid)
    eigs, vecs = model.get_eigvecs(kpts.frac)
    pjorb = np.linspace(0, model.nwan, num=model.nwan, endpoint=False, dtype=np.int) if pjorb is None else pjorb
    vecs = vecs[:, pjorb, :]
    return eigs, vecs


def proj2file(eigs, vecs, efile='eig.dat', ufile='udmft.dat', natm=1, occ=1.0, ef=0.0):

    nkpt = eigs.shape[0]
    nbnd = eigs.shape[1]
    norb = vecs.shape[1]

    # These values maybe remove in the future
    ll = 2
    nspin = 2
    nemin = 1
    nemax = nbnd
    kwt = 1.0/nkpt

    with open(efile, 'w') as fwrite:
        fmt = "{:5d}{:5d}{:5d}{:5d}{:5d}\n"
        fwrite.write(fmt.format(nkpt, natm, norb, ll, nspin))
        fmt = "{:5d}{:5d}{:5d}{:10.5f}{:10.5f}\n"
        fwrite.write(fmt.format(nemin, nemax, nbnd, occ, ef))
        fwrite.write("\n")
        for ikpt, band in enumerate(eigs):
            fwrite.write("{:5d}{:20.10f}\n".format(ikpt+1, kwt))
            for ibnd, eig in enumerate(band):
                fwrite.write("{:5d}{:20.10f}\n".format(ibnd+1, eig))
            fwrite.write("\n\n")

    with open(ufile, 'w') as fwrite:
        fmt = "{:5d}{:5d}{:5d}{:5d}{:5d}\n"
        fwrite.write(fmt.format(nkpt, natm, norb, ll, nspin))
        fmt = "{:5d}{:5d}{:5d}{:10.5f}{:10.5f}\n"
        fwrite.write(fmt.format(nemin, nemax, nbnd, occ, ef))
        fwrite.write("\n")
        fmt = "{:5d}{:5d}{:20.10f}{:20.10f}\n"
        for ikpt in range(nkpt):
            fwrite.write("{:5d}{:20.10f}\n".format(ikpt+1, kwt))
            for iorb, ibnd in itertools.product(range(norb), range(nbnd)):
                fwrite.write(fmt.format(iorb+1, ibnd+1, vecs[ikpt, iorb, ibnd].real, vecs[ikpt, iorb, ibnd].imag))
            fwrite.write("\n\n")
    return


def plotbands(kpth, eigs, file='bands.pdf', dpi=100, ylim=None, ticks=None, labels=None, linewidth=1.0,
              color='red', figsize=(6, 6), ef=0.0):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    for band in eigs:
        ax.plot(kpth, band-ef, linewidth=linewidth, color=color)
    ax.set_xlim(kpth[0], kpth[-1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    if ticks is not None:
        ax.set_xticks(ticks)
    if ticks is not None and labels is not None:
        ax.set_xticklabels(labels)
    ax.axhline(0.0, linestyle='--', color='grey')
    for tick in ticks[1:-1]:
        ax.axvline(tick, linestyle='--', color='grey')
    plt.savefig(file, dpi=dpi)
    return ax


if __name__ == "__main__":

    pjorb = [0, 1, 2, 3, 4, 5]
    eigs, vecs = gen_proj_from_string(w90file='srvo3_hr.dat', spincopy=True, pjorb=pjorb, kgrid=(20, 20, 20))
    proj2file(eigs, vecs, ef=12.5867)
