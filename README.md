![alt text](https://github.com/scarrazza/smpdf/raw/master/extra/logo.png "Logo")

Specialized Minimal PDFs tools.

##Citation

If you use this code in your publication, please cite:

*Specialized minimal PDFs for optimized LHC
calculations* [http://arxiv.org/abs/1602.00005](http://arxiv.org/abs/1602.00005) 

##Existing results

The compressed sets quoted in the paper above can be found here:

<http://pcteserver.mi.infn.it/~nnpdf/smpdf_output/>


The APPLGrids used to generate these sets can be found here:

<http://pcteserver.mi.infn.it/~nnpdf/smpdf-applgrids/>

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
executed and give valid results:

````Shell
$ lhapdf-config
$ applgrid-config
```
The default Applgrid installation **will not work**: **All** header files must 
be copied to the appromiate include path, and not only those copied by
`make install`, which is insufficient. For example, for a default
Ubuntu installation:

````Shell
applgrid-1.4.70/src $ sudo cp *.h /usr/local/include/appl_grid/
````

Python 3.4+ is also required.

In order to manage the Python dependencies,
[Anaconda](https://store.continuum.io/cshop/anaconda/) (with a Python  
3.4+ environment) is the recommended way. After setting it up,
cd into the directory of SMPDF and:

````Shell
$ conda install --file environment.yml
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
$ smpdf --help
```

for the full list of options. Note that you will probably want to
run with the `--use-db` flag.

Generating documentation
------------------------

[Sphinx](http://sphinx-doc.org/) version 1.3.1 or greater is required. 

In order to build the documentation, go into the `docs` directory, and
type:

````bash
make api
make html
````

You can use the other build targets allowed by Sphinx. Note that the
`--help` flag provides an overview of the most useful options.
