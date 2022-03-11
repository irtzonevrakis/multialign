# multialign.py: A wrapper script for the MDAnalysis structure aligner that supports multiprocessing.

# Requirements

* python >= 3.6
* [MDAnalysis](https://www.mdanalysis.org/)
* [numpy](https://numpy.org/)
* [pymp](https://github.com/classner/pymp)

Install everything with:

```
pip install MDAnalysis pymp-pypi numpy
```

# Usage

```
find /path/to/decoys/ -name '*.pdb' > pdb.lst
MKL_NUM_THREADS=1 /path/to/multialign.py --pdblist pdb.lst --outdir /path/to/aligned/decoys
```

Replace `MKL_NUM_THREADS` with `OPENBLAS_NUM_THREADS` if your numpy uses openblas instead. In any case, the script should emit a warning if the environment has not been set up correctly.

The number of threads can be customized using the `--num_threads` parameter. The default, `-1`, will use all available threads, as determined by calling `os.sched_getaffinity(0)`.

# Copyright/License

```
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
```

For more details, see the [LICENSE](LICENSE) file in this repository. For the licenses of the libraries used by this script, see the [LICENSE_COMPONENTS](LICENSE_COMPONENTS) file.
