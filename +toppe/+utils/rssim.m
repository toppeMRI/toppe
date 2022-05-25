function rss = rssim(ims, decimation)
% function rss = rssim(ims, [decimation])
%
% Recon
%
% Input
%  ims         [nx ny (nz) nCoils]
%  decimation  int     effective dwell time is decimation*4us
%
% Output 
%  rss   [nx ny (nz)]  root-sum of squares coil-combined image

sz = size(ims);

if nargin < 2
    decimation = 1;
end

nd = ndims(ims);
rss = sqrt(sum(abs(ims).^2, nd));

nx = size(rss,1);
rss = rss(end/2+((-nx/decimation/2):(nx/decimation/2-1))+1,:,:); 

