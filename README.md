![alt text](https://github.com/scarrazza/smpdf/raw/master/extra/logo.png "Logo")

Specialized Minimal PDFs tools.

## Download

Clone the master development repository by running the following command:

```Shell
$ git clone https://github.com/scarrazza/smpdf.git
```

Or download the compressed archive without using git
[here](https://github.com/scarrazza/smpdf/archive/master.zip).

## Installation

SMPDF requires [APPlgrid](https://applgrid.hepforge.org/) and
[LHAPDF](https://lhapdf.hepforge.org/). Make sure the following commands can be
executed and give valid reesults:

````Shell
$ lhapdf-config
$ applgrid-config
```

Python 3.4+ is also required.

In order to manage the Python dependencies,
[Anaconda](https://store.continuum.io/cshop/anaconda/) (with a Python  
3.4+ environment) is the recommended way. After setting it up, do:

````Shell
$ conda install conda-env
```

cd into the directory of SMPDF and:

````Shell
$ conda env update -n root -f environment.yml
```

Once all dependencies are satisfied, run:


````Shell
$ python setup.py install
```

Note that the installer script does _not_ check for dependencies. This will
install the `smpdf` program and the `smpdflib` Python library in the appropiate
paths.

## Usage

The configuration is specified in YAML files (see `examples/`). A file consists
of a list of observables (file paths to applgrids and specification of the
perturbative order, where 0 is LO and 1 NLO), a list of PDF sets (valid LHAPDF
grids in the correct path) and a list of actions (see `--help` for a
description).  

See

````Shell
$ python smpdf.py --help
```

for the full list of options. Note that you will probably want to
run with the `--use-db` flag.

