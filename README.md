# TOPPE Matlab files and user guide

This file is: https://github.com/toppeMRI/toppe/README.md

TOPPE development is primarily being done in a private Gitlab repo at the U of Michigan fMRI lab; this repo tries to stay up to date with that.

However, we want this to be a community effort and you are most welcome to propose changes! The suggested way to do so is to create a development branch and propose a merge of this branch with the master branch.


## +toppe

Matlab package (namespace) with basic toppe functions such as writemod.m, plotmod.m, playseq.m.

+toppe: core/basic functions

+toppe/+utils/: scripts for loading raw data from P-files, and other scripts.

+toppe/+utils/+rf/: make slice-selective SLR pulse with makeslr.m

+toppe/+utils/+spiral/: create stack-of-spirals readout module, and reconstruct with reconSoS.m


## examples 

Complete examples to help you get started.



## ./pulseq/

Matlab scripts for converting to/from Pulseq file format. WIP.

