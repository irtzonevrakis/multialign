#!/usr/bin/env python3
''' multialign.py: A wrapper script for the MDAnalysis structure aligner that
                   supports multiprocessing.

    Copyright (C) 2022, The Institute of Computer Science-FORTH /
                        Ioannis-Rafail Tzonevrakis

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''

from __future__ import print_function
import argparse, os, pymp
import numpy as np

from MDAnalysis import Universe
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

get_rmsd = lambda mobile, ref: rmsd(mobile.select_atoms('name CA').positions, 
                                    ref.select_atoms('name CA').positions)

def align_parallel(fnames, ref, num_threads):
    max_threads = len(os.sched_getaffinity(0))
    has_mkl = 'mkl' in ''.join(np.__config__.get_info('blas_mkl_info')\
                               ['libraries'])
    var_name = 'OPENBLAS_NUM_THREADS'
    if has_mkl:
        var_name = 'MKL_NUM_THREADS'
    suboptimal_np_threads = False
    if var_name in os.environ.keys():
        if os.environ[var_name] != '1':
            suboptimal_np_threads = True
    else:
        suboptimal_np_threads = True
    if suboptimal_np_threads:
        print(f'WARN: Running with {var_name} != 1. Sub-optimal performance '+
              f'may ensue. We recommend re-running the script with {var_name}'+
              ' set to 1.')
    if num_threads > max_threads:
        print(f'WARN: Number of threads {num_threads} larger than maximum '+
              f'alotted threads {max_threads}. Setting to {max_threads} instead.')
        num_threads = max_threads
    if num_threads == -1:
        num_threads = max_threads
    print(f'Running with {num_threads} threads.')
    fnames = np.array(fnames)
    shmem = pymp.shared.array((len(fnames),), dtype=fnames.dtype)
    shmem = fnames
    with pymp.Parallel(num_threads) as p:
        for i in p.range(1, len(fnames)):
            mobile = Universe(shmem[i])
            print(f'Decoy {fnames[i]}, RMSD to ref before alignment: '+
                  f'{get_rmsd(mobile, ref)}')
            align.AlignTraj(mobile, ref, in_memory=True).run()
            mobile.select_atoms('all').write(os.path.join(outdir, 
                                             os.path.basename(shmem[i])))
            print(f'Decoy {fnames[i]}, RMSD to ref after alignment: '+
                  f'{get_rmsd(mobile, ref)}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--pdblist', type=str, required=True,
                        help='Path to list of PDBs. The first PDB is used as '+
                             'reference structure.')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Path to directory to store aligned PDBs in.')
    parser.add_argument('--num_threads', type=int, required=False, default=-1,
                        help='Number of threads to use (defaults to maximum '+
                              'available).')
    args = parser.parse_args()
    pdblist = args.pdblist
    outdir = args.outdir
    num_threads = args.num_threads
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fnames = []
    with open(pdblist) as fp:
        for fname in fp:
            fnames.append(fname.strip())

    print(f'Using {fnames[0]} as reference')
    ref = Universe(fnames[0])
    ref.select_atoms('all').write(os.path.join(outdir, 
                                  os.path.basename(fnames[0])))

    align_parallel(fnames, ref, num_threads)
