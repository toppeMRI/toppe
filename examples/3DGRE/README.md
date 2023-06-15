# 3D SPGR sequence example

Design a 3D RF-spoiled GRE sequence (FLASH/SPGR/T1-FFE) 
in TOPPE format for execution on GE scanners.

## Usage

1. Create scan files (gre.tar)
    ```
    >> gre;
    ```

2. Copy the contents of gre.tar to /usr/g/research/pulseq/v5/gre/
   on the scanner host computer.
    ```
    $ pwd
    /usr/g/research/pulseq/v5/gre/
    $ tar xf gre.tar
    ```

3. Move toppe0.entry to /usr/g/research/pulseq/v5/:
    ```
    $ pwd
    /usr/g/research/pulseq/v5/gre/
    $ mv toppe0.entry ..
    ```

4. Prescribe the TOPPE interpreter with the following settings:
    ```
    opuser1 = 0   (selects toppe0.entry as the entry point)
    oprbw = 31.25 (such that ADC/dwell time = 16us).
    ```
    Other settings are as shown in the
    [TOPPE EPIC source code repository](https://github.com/jfnielsen/TOPPEpsdSourceCode/).

5. (optional) Reconstruct and display:
    ```
    recongre;
    ```



## Troubleshooting

* See the
[TOPPE EPIC source code repository](https://github.com/jfnielsen/TOPPEpsdSourceCode)
for common failure modes and troubleshooting tips.

* If using `toppe.utils.rf.makeslr()` to create tipdown.mod 
and you're not able to compile the mex files for the SLR toolbox, 
try using the Matlab scripts instead by copying them to somewhere in your path:
    ```
    $ cp ../../+toppe/+utils/+rf/+jpauly/matlab/ .
    ```
