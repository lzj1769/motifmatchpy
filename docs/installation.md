# Installation

## Install a release

Python **3.12 or newer** is required.

```sh
python -m pip install motifmatchpy
motifmatchpy --version
motifmatchpy --help
```

With uv, add it to your project:

```sh
uv add motifmatchpy
uv run motifmatchpy --version
```

The package depends on NumPy and Numba. It does not build or load the MOODS C++
core at runtime. These dependencies contain native components; installation
without a compiler depends on compatible wheels being available for your Python
version and platform. The first scan may take longer while Numba compiles its
kernels; subsequent calls reuse them.

## Install the current source

```sh
git clone https://github.com/lzj1769/motifmatchpy.git
cd motifmatchpy
python -m pip install .
```

This also provides the `docs/examples/` files used throughout the guide.
Installing the main package does not install the development-only C++ reference.
For editable development and tests, see [contributing](development.md).

## Check the installation

```sh
python -c "import motifmatchpy as mm; print(mm.__version__, mm.MOODS_VERSION)"
python -m motifmatchpy.cli --help
```

If `motifmatchpy` is not on your shell's PATH, activate the environment containing
the installation, or use `python -m motifmatchpy.cli` in place of the executable.
