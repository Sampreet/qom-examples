# Examples for the Quantum Optomechanics Toolbox

[![Toolbox](https://img.shields.io/badge/qom-1.0.0-red?style=for-the-badge)](https://sampreet.github.io/qom-docs)
[![Last Commit](https://img.shields.io/github/last-commit/sampreet/qom-examples?style=for-the-badge)](#)

> A set of notebooks and scripts to demonstrate the usage of the quantum optomechanics toolbox!

## Notebooks

* [Dynamics of an End-mirror Optomechanical System](notebooks/systems/em_00.ipynb)

## Structure of the Repository

```
ROOT_DIR/
|
├───data/
│   ├───bar/
│   │   ├───foo_baz_xyz.npz
│   │   └───...
│   └───...
|
├───notebooks/
│   ├───bar/
│   │   ├───foo.ipynb
│   │   └───...
│   └───...
|
│───scripts/
│   ├───bar/
│   │   ├───foo_baz.py
│   │   └───...
│   └───...
|
├───systems/
│   ├───__init__.py
│   ├───Foo.py
│   └───...
│
├───.gitignore
├───CHANGELOG.md
└───README.md
```

## Execution

### Installing Dependencies

The project requires `Python 3.8+` installed via the [Anaconda distribution](https://www.anaconda.com/products/individual). 
An extensive guide to set up your python environment same can be found [here](https://sampreet.github.io/python-for-physicists/modules/m01-getting-started/m01t01-setting-up-python.html).

Once the installation is complete and `conda` is configured, it is preferable to create a new conda environment (say `qom`) and activate it using:

```bash
conda create -n qom python=3
conda activate qom
```

This project uses [The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom) which can be installed via Python Package Index using `pip` by executing:

```bash
pip install -i https://test.pypi.org/simple/ qom
```

Alternatively, [clone the repository](https://github.com/Sampreet/qom) or [download the sources](https://github.com/Sampreet/qom/archive/refs/heads/master.zip) as `.zip` and extract the contents.
Now, execute the following from *outside* the top-level directory, `ROOT_DIR`, inside which `setup.py` is located:

```bash
pip install -e ROOT_DIR
```

### Running the Scripts

To run the scripts, navigate *inside* the top-level directory, `ROOT_DIR`, and execute:

```bash
python scripts/bar/baz.py
```

Here, `bar` is the name of the folder inside `scripts` and `baz.py` is the name of the script (refer to the repository structure).