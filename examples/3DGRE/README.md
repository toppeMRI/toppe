## 3D SPGR sequence example

Complete example, including:

+ Create .mod files using the TOPPE Matlab library
+ Create scanloop.txt, and modules.txt
+ Display the sequence in movie (loop) mode
+ Create seqstamp.txt, needed to execute the sequence on scanner

Usage:
```
>> gre;
```

Issue:
If you're not able to compile the mex files for the SLR toolbox, try using the Matlab scripts instead by copying them to somewhere in your path:
```
$ cp ../../+toppe/+utils/+rf/+jpauly/matlab/ .
```
