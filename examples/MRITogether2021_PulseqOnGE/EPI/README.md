# Variable flip angle 3D multi-shot EPI SPGR/FLASH sequence

Create the sequence files:
```
>> epi;
```

Display sequence in loop/movie mode:
```
>> toppe.playseq(4, seq.sys, 'tpause', 0.05);
```

Load P-file and display images:
```
>> loadEPI('P,epi.7', 'readout.mod');
```
