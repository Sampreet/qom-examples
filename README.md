# Examples for the Quantum Optomechanics Toolbox

[![Toolbox Version](https://img.shields.io/badge/qom->1.0.1-red?style=for-the-badge)](https://sampreet.github.io/qom-docs)
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

[The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom) requires `Python 3.8+` installed, preferably via the [Anaconda distribution](https://www.anaconda.com/download).
Once `Anaconda` is set up, create and activate a new `conda` environment using:

```bash
conda create -n qom python
conda activate qom
```

The toolbox relies primarily on `numpy` (for fast numerical algebra), `scipy` (for numerical methods), `sympy` (for symbolic algebra), `seaborn` (for color palettes) and `matplotlib` (for plotting results).
These libraries can be installed using:

```bash
conda install matplotlib numpy scipy sympy seaborn
```

***Note: To run the GUI modules, `pyqt` should be installed separately.***

Once the dependencies are installed, download the [repository of the toolbox](https://github.com/Sampreet/qom) as `.zip` and extract the contents.
Now, execute the following from *outside* the top-level directory, `ROOT_DIR`, inside which `setup.py` is located (refer to the file structure of the repository [here](https://github.com/sampreet/qom/blob/master/CONTRIBUTING.md)):

```bash
pip install -e ROOT_DIR
```

The corresponding documentation is available [here](https://sampreet.github.io/qom-docs).

### Running the Scripts

To run the scripts, navigate *inside* the top-level directory of this repository, `ROOT_DIR`, and execute:

```bash
python scripts/bar/foo_baz.py
```

Here, `bar` is the name of the folder inside `scripts` and `foo_baz.py` is the name of the script (refer to the repository structure).