#!/usr/bin/env python
# coding=UTF-8
'''
@Author: Jianzhou Zhao
@Date: 2019-07-15 20:25:49
@LastEditors: Jianzhou Zhao
@LastEditTime: 2019-07-16 11:12:36
@Description: This code is used for run DMFT loop by combining KSUM and IQIST.
'''


def check_file_status(efile, ufile, mfile, sfile):
    from os.path import isfile
    if not isfile(efile):
        exit('Eigenvalue file is not exist.')
    if not isfile(ufile):
        exit('Projection file is not exist.')
    if not isfile(mfile):
        exit('DMFT main input file is not exist.')
    if not isfile(sfile):
        exit('Solver main input file is not exist.')
    return


def mainloop(maxloop=1, wrap_root=None, ksumcmd=None, solvercmd=None):

    import subprocess

    status = 0

    for iloop in range(maxloop):

        ksum_output = 'ksum_'+str(iloop)+'.out'
        solver_output = 'ctqmc_'+str(iloop)+'.out'

        # KSUM run
        print(iloop+1, ksumcmd)
        process = subprocess.Popen(ksumcmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = process.communicate()
        with open(ksum_output, 'w') as output:
            output.write(str(out, encoding="utf8"))
        if err is not None:
            statu = 1
            break

        # Solver run
        print(iloop+1, solvercmd)
        process = subprocess.Popen(solvercmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = process.communicate()
        with open(solver_output, 'w') as output:
            output.write(str(out, encoding="utf8"))
        if err is not None:
            status = 2
            break

    return status


if __name__ == "__main__":

    import os

    efile = 'eig.dat'
    ufile = 'udmft.dat'
    mfile = 'dmft.main.in'
    sfile = 'solver.ctqmc.in'
    maxloop = 5

    wrap_root = os.path.dirname(os.path.realpath(__file__))
    mpicmd = 'mpirun'
    mpiopt = '-n 4'
    ksumcmd = mpicmd + ' ' + mpiopt + ' ' + wrap_root+'/../ksum/ksum'
    solvercmd = mpicmd + ' ' + mpiopt + ' ' + wrap_root+'/../solver/gardenia/ctqmc'

    check_file_status(efile, ufile, mfile, sfile)

    status = mainloop(maxloop=maxloop, wrap_root=wrap_root, ksumcmd=ksumcmd, solvercmd=solvercmd)
