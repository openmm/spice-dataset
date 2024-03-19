# Dataset Submission

This directory contains scripts for submitting datasets to be calculated on QCFractal.  To use them, you must have
QCPortal 0.5 or later installed.  The input is an HDF5 file describing the samples in the dataset to create.  These
files are found in the other directories in this repository.

`submit.py` submits a new dataset.  The first argument is the name of the dataset to create.  It is followed by one or
more HDF5 files defining the samples.

`checkStatus.py` reports the status of performing the calculations for a dataset.  It takes a single argument, the name
of the dataset.

`resetErrors.py` resets the status of calculations that have failed so they will be attempted again.  It takes a single
argument, the name of a dataset.  Any record whose status is "error" is reset to "waiting".