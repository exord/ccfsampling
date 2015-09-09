# ccfsampling
A test of the impact of CCF sampling in RV diagnostics precision.

Prerequisites:
- PASTIS_v6
- numpy
- pylab

The main file in the repository is test_ccfsampling.py.

It can be run directly as a python script, and produces plots of the relative and absolute precision of the Bisector, vspan and BiGauss diagnostics as a function of the binning level of the CCF. It also produces a plot of the time requiered for each computation.

Usage
-----
Directly as python script from the shell.

$ python test_ccfsampling.py

Or from within python itself.

\>\>\> execfile('test_ccfsampling.py')

