# TOPPE Matlab toolbox

TOPPE is a modular framework for rapid prototyping of pulse sequences on General Electric MRI scanners.
For an overview, see <https://toppemri.github.io/>.

TOPPE can be used on its own, or as the GE interpreter for [Pulseq](http://pulseq.github.io/); 
the latter allows MR sequences to be ported across vendor platforms (currently, Siemens and GE are supported).

This repository contains various basic functions for writing and reading TOPPE scan files, that are
then simply copied to the scanner host computer and loaded by the TOPPE interpreter for execution.
The interpreter (EPIC source code) is hosted in a 
[separate repository](https://github.com/jfnielsen/TOPPEpsdSourceCode) -- for access, contact jfnielse@umich.edu.


## Example sequence: 3D GRE

The [examples/3DGRE](examples/3DGRE) folder contains the Matlab function `gre.m`, that
creates a 3D RF-spoiled gradient-echo sequence:

1. Get this toolbox
```
$ git clone git@github.com:toppeMRI/toppe.git
```
1. Create the GRE scan files.
In Matlab:
```
>> cd examples/3DGRE
>> addpath ../..   % folder containing the +toppe package
>> gre;
```
This will create the files `modules.txt`, `scanloop.txt`, `seqstamp.txt`, and two `.mod` files (`tipdown.mod` and `readout.mod`).

1. Place those files in the folder **/usr/g/research/pulseq/gre/** on scanner host.
1. Create a text file `toppe<CV1>.meta` in **/usr/g/research/pulseq/** with the following contents:
```
/usr/g/research/pulseq/gre/
modules.txt
scanloop.txt
tipdown.mod
readout.mod
seqstamp.txt
```
(This file is also provided in the examples/3DGRE folder, with the name toppe0.meta.)
Here, `CV1` is the non-negative integer you enter in the scan prescription (see below).
Note the **trailing backslash** in the path on the first line.

The folder /usr/g/research/pulseq/ should now look something like this:
![Example directory contents](resource/folder.png)

For a more detailed explanation of these files, see [Files.md](Files.md).


### Running the tv4 interpreter

1. On the Imaging Options screen, prescribe an axial 2D GRE scan and type `/usr/g/research/pulseq/tv4` in the PSD name box:  
   ![Imaging options screen](resource/imagingOptions_DV26.png)
1. On the main scan prescription page, use the following settings: 
   1. Slice Thickness = 1.0
   1. Spacing = 0.0
   1. Set S/I Start and End locations to `S120` and `I120`. This will ensure that 
   the number of slices _is greater than maxslice in scanloop.txt_, which is a requirement 
   (but the exact number of slices is not important, and you may want to use a large number so the same protocol can be used for other TOPPE scans!)  
   1. Freq. Dir = L/R
   1. Intensity correction: NONE
   1. Intensity Filter: NONE
   1. Shim: whatever setting you like
   ![RxScreen](resource/RxScreen1_DV26.png)
1. In the `Advanced` tab, use the following settings:
   1. CV1 = Non-negative integer that determines the `.meta` file that serves as the entry point for TOPPE.
   This value is used to select between different TOPPE scans, and is saved with the protocol.
   The default is 0, in which case the scan entry point is the file `/usr/g/research/pulseq/toppe0.meta`.
   1. CV2 = 1 to save a Pfile, 0 to only save the ScanArchive file.
   1. CV3 = 64. This value can be changed to optimize sequence timing (advanced use) (us).
   1. CV4 = 200. Sequence ssi time (us).
   ![RxScreen](resource/RxScreen2_DV26.png)
1. Press `Save Rx`
1. From the pull-down menu near the Scan button, select `Download` (this will take a few seconds), then `Auto Prescan`  
   ![Download](resource/download.png)
1. Scan. The scan clock should show the scan time listed in scanloop.txt.
1. For interactive slice planning, one option is to use the Java GUI in https://github.com/toppeMRI/SlicePlanner.


## Troubleshooting

Under construction


## Version history

I have been sloppy with Git tags/versions so here is a brief history.
The current version is toppev3 which we have used for over a year now and which seems stable.

toppev2b: like toppev2, but should guarantee correct b1 scaling

toppev2c: like toppev2b, but shows scan time on console (see entry in scanloop.txt)

toppev2d: like toppev2c, except reads a file 'toppescan<N>.meta', where <N> = value of usercv,
          that lists what files to use. Also lists the "primary" RF and readout modules, that are used to
          populate certain internal EPIC fields related to b1 scaling, SAR, acquisition filter, etc.   

toppev3: like toppev2d, except scanloop.txt contains additional columns that specify 3D rotation matrix.

tv4: Adds system and patient safety checks. See ./tv4/README.md.
