# Dataset Submission

This directory contains scripts for submitting datasets to be calculated on QCFractal.  To use them, you must have
QCPortal 0.5 or later installed.  The input is an HDF5 file describing the samples in the dataset to create.  These
files are found in the other directories in this repository.

`submit.py` takes two arguments: the name of the dataset to create, and the HDF5 file defining the samples.

`checkStatus.py` reports the status of performing the calculations for a dataset.  It takes a single argument, the name
of the  dataset.