function C = constants

C.NULL = 0;
C.TRAP = 1;
C.ARBITRARY = 2;
C.ADC = 1;
C.MAXIAMP = 2^15-2;  

% internal constants used to avoid writing floats to .mod file
C.RFSCALE = 1e4;  
C.GSCALE = 3000;  
