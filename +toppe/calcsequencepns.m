function p = calcsequencepns(sysGE, varargin)

% defaults
arg.timeRange = [0 inf];  % the whole sequence
arg.plot = true;
arg.print = true;

arg = toppe.utils.vararg_pair(arg, varargin);

[~, gx, gy, gz] = toppe.plotseq(sysGE, 'timeRange', arg.timeRange, 'doDisplay', false);

raster = sysGE.raster*1e-6;   % s
grad(1,:) = gx(:)'*1d-2;    % T/m
grad(2,:) = gy(:)'*1d-2;
grad(3,:) = gz(:)'*1d-2;
p = toppe.pns(grad, sysGE.gradient, 'gdt', raster, 'plt', arg.plot, 'print', arg.print);

