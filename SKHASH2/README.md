# SKHASH

Python package for earthquake focal mechanism inversions.

Author: Robert Skoumal, U.S. Geological Survey | rskoumal@usgs.gov

This project contains Python code to compute focal mechanism solutions using first-motion polarities (traditional, consensus, and/or imputed) and S/P ratios (traditional and/or consensus).

## How to Use
1. Clone the project
```
git clone https://code.usgs.gov/esc/SKHASH.git
```
2. Python 3.8+ versions are supported. Ensure the [NumPy](https://numpy.org) and [pandas](https://pandas.pydata.org) packages are available in your virtual environment. If you wish to create the automatic beachball plots, [Matplotlib](https://matplotlib.org) will also be needed.
```
pip3 install numpy pandas
pip3 install matplotlib # Only needed if creating beachballs

# Alternatively:
pip3 install -r requirements.txt
```
3. Run SKHASH and provide the path to your desired control file.
```
python3 SKHASH.py examples/hash1/control_file.txt
```
4. Refer to the [manual](https://code.usgs.gov/esc/SKHASH/-/blob/main/SKHASH_manual.pdf) for additional information about running the code.

## Fortran subroutine (completely optional)
By default, SKHASH will compute mechanisms using the Python routine. However, to speed up the grid search, you can choose to take advantage of an included Fortran subroutine. If the Python C/API wrapper does not already exist, SKHASH will automatically create the wrapper with the user's permission when the Fortran subroutine is used.

To use this Fortran routine, add the following lines to your control file:
```
$use_fortran
True
```

A fortran compiler will be needed on the user's machine. If one does not exist, here are some examples on how to get one:
```
# macOS Homebrew example
brew install gcc

# Ubuntu example
apt install gfortran
```

If you're using Python 3.12+, the [meson](https://pypi.org/project/meson-python) and [ninja](https://pypi.org/project/ninja) packages are also needed.
```
pip3 install meson ninja
```

Note that if you are using macOS and receive an error, you may be missing the Command Line Tools package. To install:
```
xcode-select --install
```

## Citation
Please cite our paper if you use anything in this project:

- Skoumal, R.J., Hardebeck, J.L., Shearer, P.M. (_in press_). SKHASH: A Python package for computing earthquake focal mechanisms. _Seismological Research Letters_.

Significant portions of this algorithm are based on [HASH](https://www.usgs.gov/node/279393):

- Hardebeck, J.L., & Shearer, P.M. (2002). A new method for determining first-motion focal mechanisms. _Bulletin of the Seismological Society of America_, 92(6), 2264-2276. https://doi.org/10.1785/0120010200

- Hardebeck, J.L., & Shearer, P.M. (2003). Using S/P amplitude ratios to constrain the focal mechanisms of small earthquakes. _Bulletin of the Seismological Society of America_, 93(6), 2434-2444. https://doi.org/10.1785/0120020236

## License and Disclaimer
[License](https://code.usgs.gov/esc/SKHASH/-/blob/main/LICENSE.md): This project is in the public domain.

[Disclaimer](https://code.usgs.gov/esc/SKHASH/-/blob/main/DISCLAIMER.md): This software is preliminary or provisional and is subject to revision.
