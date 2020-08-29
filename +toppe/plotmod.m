function plotmod(fname, varargin)
% Plot TOPPE .mod file.
%
% function plotmod(fname)
%
% Inputs:
%   fname      [string]    .mod file name, or 'all' (print all .mod files in current folder)
% Options:
%  plotPNS     [boolean]   default: true
%  printPNS    [boolean]   print pns numbers to console. Default: false.
%  gradcoil    [string]    gradient subsystem (see pns.m). Default: 'xrm'

arg.plotPNS = true;
arg.printPNS = false;
arg.gradcoil = 'xrm';

arg = toppe.utils.vararg_pair(arg, varargin);

if strcmp(fname, 'all')
	f = dir('*.mod');
	for ii = 1:length(f)
		sub_plotmod(f(ii).name, arg);
	end
else
	sub_plotmod(fname, arg);
end

return;

function sub_plotmod(fname, arg)

[b1,gx,gy,gz] = toppe.readmod(fname);
rho = abs(b1);
theta = angle(b1);

figure;

% PNS
if arg.plotPNS
	gdt = 4d-6;           % raster time (sec)
	nwavs = size(gx,2);   % number of waveforms in .mod file
	for ii = 1:nwavs
		clear grad;
		grad(1,:) = gx(:,ii)'*1d-2;    % T/m
		grad(2,:) = gy(:,ii)'*1d-2;
		grad(3,:) = gz(:,ii)'*1d-2;
		plt = false;           % plot output
		% [p.PThresh, p.pt, p.PTmax, p.gmax, p.smax] = toppe.pns(grad, arg.gradcoil, 'gdt', gdt, 'plt', false, 'print', arg.printPNS);
		[p.PThresh] = toppe.pns(grad, arg.gradcoil, 'gdt', gdt, 'plt', false, 'print', arg.printPNS);
		pns.val(:,ii) = p.PThresh(:);
	end

	subplot(4,3,11);
	%t = (0:(size(p.pt,2)-1))*gdt*1e3;     % time (ms)
	t = (0:(size(pns.val,1)-1))*gdt*1e3;     % time (ms)
	%plot(t,p.pt(:,:,1),'',t,p.PThresh(:,:,1),'r--',...
   %     [t(1) t(end)],[100 100],'m:',[t(1) t(end)],[80 80],'m:',...
   %     [t(1) t(end)],-[100 100],'m:',[t(1) t(end)],-[80 80],'m:');
	plot(t, pns.val, '', ...
		[t(1) t(end)],[100 100],'m:', ...
		[t(1) t(end)],[80 80],'m:');
	xlabel('time [ms]'); ylabel('PNS [% of threshold]');
	tmp = 1.05*max([p.PThresh(:,:,1) 100]);
	grid on; % axis([0 t(end) 0 tmp]);
end

nt = size(b1,1);
dt = 4e-3;  % ms
T = linspace(dt/2,nt*dt-dt/2,nt);

% rf
subplot(437); sub_plot(T, rho);   ylabel('abs(rf) G');
xlabel('time (msec)');
subplot(438); sub_plot(T, theta); ylabel('angle(rf) rad');

% gradient waveform
subplot(431); sub_plot(T,gx);    ylabel('gx G/cm');   
try
	sgtitle(fname);
catch
	title(fname);
end
subplot(432); sub_plot(T,gy);    ylabel('gy G/cm');
subplot(433); sub_plot(T,gz);    ylabel('gz G/cm');

% gradient slew
subplot(434); sub_plot(T(2:end), diff(gx)/dt);  ylabel('x slew  G/cm/ms');
subplot(435); sub_plot(T(2:end), diff(gy)/dt);  ylabel('y slew  G/cm/ms');
subplot(436); sub_plot(T(2:end), diff(gz)/dt);  ylabel('z slew  G/cm/ms');
g = [gx gy gz];
slew = sqrt(sum((diff(g,1)/dt).^2,2));
subplot(4,3,10); sub_plot(T(2:end), slew, 'Color', 'r');  ylabel('combined slew  G/cm/ms');

return;

function sub_plot(T,wav,varargin)

arg.Color = 'b';
arg.LineWidth = 1.0;
arg = toppe.utils.vararg_pair(arg, varargin);

h = line(T,wav);
set(h, 'LineWidth', arg.LineWidth);
set(h, 'Color', arg.Color);

return;

