function C = constants

C.NULL = 0;
C.TRAP = 1;
C.ARBITRARY = 2;
C.ADC = 1;
C.max_pg_iamp = 2^15-2;  

% internal constants used to avoid writing floats to .mod file
C.rfscale = 1e4;  
C.gscale = 3000;  
