% perform centered fft/ifft
%
% function K = cfftn(D, direction)

function K = cfftn(D, direction)

if nargin<2
	error('toppe.utils.fftn: provide direction as either ''forward'' or ''inverse'' ');
end

if strcmp('forward', direction)
	K = fftshift(fftn(fftshift(D)));
else
	K = fftshift(ifftn(fftshift(D)));
end
