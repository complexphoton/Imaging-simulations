# Imaging simulations

## Introduction
This repository contains the code for performing large-scale full-wave simulations of five optical imaging methods by the scattering matrix: 
- Scattering matrix tomography ([SMT](https://arxiv.org/abs/2306.08793))
- Reflectance confocal microscopy (RCM)
- Optical coherence tomography (OCT)
- Optical coherence microscopy (OCM)
- Interferometric synthetic aperture microscopy (ISAM)

In principle, our approach can model any scattering-based imaging method. Such numerical modeling can be helpful for developing new imaging methods by providing the ground truth, the flexibility to tailor the system and the imaging scheme, and the ease of comparing different methods.

## Installation

### Prerequisites
In order to run the code, you will need to install the following dependencies.

- System setup: [MESTI](https://github.com/complexphoton/MESTI.m), [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) (optional, a package for producing the matrix ordering)
- Reflection matrix computation: [MESTI](https://github.com/complexphoton/MESTI.m), [MUMPS](https://mumps-solver.org/index.php)
- Image reconstruction: [MESTI](https://github.com/complexphoton/MESTI.m), [FINUFFT](https://finufft.readthedocs.io/en/latest/)

Note that if you plan to use only one component of the code, there is no need to install dependencies of other components.

### Install Dependencies

#### MESTI

The installation of MESTI is straightforward. Simply download the source code from [here](https://github.com/complexphoton/MESTI.m/tree/main/src) and add it to the MATLAB search path.

#### METIS

Download the source code from [here](http://glaros.dtc.umn.edu/gkhome/metis/metis/download). No need to use the OpenMP version.

You can compile METIS using the default option:
```
$ make config
$ make
```

After compilation, you can install METIS by running:
```
$ make install
```
The default installation path is /usr/local, which is in the MATLAB search path.

#### MUMPS

The detailed instructions for compiling MUMPS can be found [here](https://github.com/complexphoton/MESTI.m/tree/main/mumps). Be cautious that the compilation of MUMPS and its MATLAB interface can take significant effort. 

#### FINUFFT

Download the source code from [here](https://github.com/flatironinstitute/finufft). Currently, there are two routes to compile: the new CMake route or the old GNU makefile route. 
<!--
Currently, there are two routes to compile: CMake or GNU makefile. My compilation was done using the GNU makefile route. 
-->

There is a caveat in compiling the MATLAB interface on Mac. You may get a warning of ```license has not been accepted``` from Xcode and a following error of ```no supported compiler was found```. The same error can happen to the MUMPS MATLAB interface. The simple [solution](https://finufft.readthedocs.io/en/latest/install.html#the-clang-route-default) is typing 

```
$ /usr/libexec/PlistBuddy -c 'Add :IDEXcodeVersionForAgreedToGMLicense string 10.0' ~/Library/Preferences/com.apple.dt.Xcode.plist
```
in the command line. 

In the future, we may add an option to use the [nufftn](https://www.mathworks.com/help/matlab/ref/double.nufftn.html) from MATLAB, which does not require the compilation of FINUFFT but is slower.

After installing the required dependencies above, you can simply download the code and add all folders to the MATLAB search path. Make sure that the working directory is ```/path/to/Imaging-simulations``` when you run the code.

## Getting Started
To get started, we suggest running the image reconstruction script for the system described in our paper, after downloading the precomputed reflection matrix and the corresponding ```system_data.mat``` file. Each image can be reconstructed in one minute on a MacBook Air with Apple M1 chip. Running ```recon_all.m``` will reconstruct all images.


## Usage

### System Setup

### Reflection Matrix Computation

### Image Reconstruction

### Post-processing

### Plotting

## Results
