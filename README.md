# TOPPE Matlab toolbox

This repository contains various basic functions for writing and reading TOPPE scan files, that are
then simply copied to the scanner host computer and loaded by the TOPPE interpreter for execution.
For a detailed explanation of the files that make up a TOPPE scan,
see [Files.md](Files.md).

**UPDATE Aug 2023:**
The recommended way to use TOPPE is to first design the sequence in 
[Pulseq](http://pulseq.github.io/), and then convert the .seq file to a .tar
file that can be run on a GE scanner with the TOPPE interpreter.
The interpreter (EPIC source code) is hosted in a 
[separate repository](https://github.com/jfnielsen/TOPPEpsdSourceCode) -- for access, contact jfnielse@umich.edu.
That repository also hosts the new **Pulseq on GE User Guide**.

For this reason, we have moved some of the information in this repository
(including the examples) to the 'attic' folder.
