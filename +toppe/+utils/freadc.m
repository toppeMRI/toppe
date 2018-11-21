function C = freadc(fid, len)
%freadc - read and deblank character strings
%
%  C = freadc(fid, len)
%    fid - open file handle to read
%    len - number of characters to read
%    C - string read

% Copyright (c) 2012 by General Electric Company. All rights reserved.

% Read and deblank character strings
C = deblank( char( fread( fid, [1, len], 'uchar')));
