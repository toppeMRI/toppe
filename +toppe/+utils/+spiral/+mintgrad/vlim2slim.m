%
%	Function [gmax,smax] = vlim2slim(Imax,Vmax,Rcoil,Lcoil,eta)
%
%	Function takes parameters form a voltage-limited model and
%	returns the gmax and smax of the most conservative 
%	constant-slew-rate-limited model that will satisfy the voltage
%	model.
%
%	INPUT:
%		Imax  = maximum amplifier current (A)
%		Vmax  = maximum amplifier voltage (V)
%		Rcoil = coil resistance (ohms)
%		Lcoil = coil inductance (H)
%		eta   = coil efficiency (G/cm/A)
%
%	OUTPUT:
%		gmax  = maximum gradient (G/cm)
%		smax  = maximum slew (G/cm/s)
%
%
%	B. Hargreaves, April 2003.


function [gmax,smax] = vlim2slim(Imax,Vmax,Rcoil,Lcoil,eta)

if (nargin < 5)
	error('Not enough arguments');
end;


gmax = eta.*Imax;			% Max grad is just proportional to Imax.
smax = eta.*(Vmax-Rcoil.*Imax)./Lcoil;	% Max slew assumes worst-case grad.


