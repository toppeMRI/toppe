# 2D GRE and EPI sequence examples for ISMRM 2019 software demo, Sunday 12 May.

Use with psd version 'toppev2c'. See https://toppemri.github.io/ for more info.


## 2D GRE

To generate TOPPE files:
```
gre_live_demo
```

To recon the Pfile 'Pgre.7' which was acquired with 'readout.mod':
```
gre_recon('Pgre.7', 'readout.mod');
```


## 2D EPI

To generate TOPPE files:
```
epi_live_demo
```

To recon the Pfile 'Pepi.7' which was acquired with 'readout.mod':
```
gradDelay = 2e-6;                      % seconds
oddEvenPhaseOffset = -0.2;             % radians
epi_recon('Pepi.7', 'readout.mod', gradDelay, oddEvenPhaseOffset);
```

