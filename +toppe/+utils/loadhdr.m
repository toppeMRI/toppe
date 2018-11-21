function [rdb_hdr] = loadhdr(pfile)
% function [rdb_hdr] = loadhdr(pfile)
%
% Calls read_rdb_hdr on pfile.
%
% $Id: loadhdr.m,v 1.2 2018/10/25 13:14:40 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/loadhdr.m,v $

import toppe.utils.*

fid = fopen(pfile,'r','l');
ver = fread(fid,1,'float32');
str = num2str(ver);
fprintf('Pfile version is %s\n', str);
rdbm_rev = str2double(str);
fseek(fid,0,'bof');                 % NB!
rdb_hdr = read_rdb_hdr(fid,rdbm_rev);

return;
