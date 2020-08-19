function wav = sub_trap2shape(grad, gamma)
% Convert trapezoid to (positive) shape on 4us sample boundary
% Increase rise/fall times by one sample to keep within slew limit (due to 4us sampling)
%
% Example:
%  seq = mr.Sequence();
%  seq.read('myseq.seq');
%  trap = seq.getBlock(2).gz
%  gamma = toppe.systemspecs().gamma;
%  wav = sub_trap2shape(trap, gamma);

dt = 4e-6;   % sec
gamp = abs(grad.amplitude)/100/gamma;                % Gauss/cm
wav = [ linspace(0, gamp, ceil(grad.riseTime/dt)+1)  ...
		gamp*ones(1, floor(grad.flatTime/dt)) ... 
		linspace(gamp, 0, ceil(grad.fallTime/dt)+1) ]';

return
