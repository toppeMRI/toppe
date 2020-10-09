function M = blochsim(Mi, beff, T1, T2, dt, nstep)
%
%Block equation simulator used in Doug Noll's BME516 Medical Imaging Systems Course
%
%function M = blochsim(Mi, beff, T1, T2, dt, nstep)
%
%need to supply the follwing parameters:
%
%Mi = [Mxi,Myi,Mzi]  initial magnetization in the x,y, and z directions
%
%beff = nstep x 3 array.  the columns are Beff (in Tesla) along the x,y,and z directions
%
%T1,T2 are the relaxation constants (in ms)
%
%dt step size (in ms) of simulation
%
%nstep = # of steps in simulation

gambar = 42.57e3;        % gamma/2pi in kHz/T
gam = gambar*2*pi;
T1 = dt ./ T1;		% trick - now it is the inverse!
T2 = (1 - dt ./ T2);		% trick: now it is the loss per step

%
% Put Beff into appropriate units
%
beff = beff.*(dt*gam);
M = zeros(nstep,3);
M(1,:) = Mi;
runge = 1;  % alternate method is faster and less sensitive to dt ...

if runge == 0
    % 4th order runge kutta
    m = Mi;
    for lp = 2:nstep
        k1 = cross(m,beff(lp-1,:));
        k2 = cross(m+k1./2,beff(lp-1,:));
        k3 = cross(m+k2./2,beff(lp-1,:));
        k4 = cross(m+k3,beff(lp,:));
        m = m + (k1 + 2.*k2 + 2.*k3 + k4)./6;
        m = m.*[T2 T2 (1-T1)] + [0 0 T1]; % relaxation
        M(lp,:) = m;
    end
    
else
    
    % alternate form where rotations are explicitly calculated and carried out
    % on the magnetization vector - this is a more stable solution
    for lp = 2:nstep
        B = beff(lp-1,:);
        %
        %	Compute sines & cosines of field angles:
        %	Theta = angle w.r.t positive z axis
        %	Phi   = angle w.r.t positive x axis
        %	Psi   = angle w.r.t transformed positive x axis
        %
        Bmag = sqrt(sum(B.^2));		% Magnitude of applied field
        Btrans = sqrt(B(1).^2 + B(2).^2);	% Magnitude of transverse applied field
        
        ct = 1;
        if Bmag ~= 0;
            ct = B(3) ./ Bmag;	% cos(theta)
        end
        st = sqrt(1 - ct.^2);				% sin(theta) > 0
        
        cphi = 1;
        if Btrans ~= 0;
            cphi = B(1) ./ Btrans;	% cos(phi)
        end
        sphi = sqrt(1 - cphi.^2) .* sign(B(2));	% sin(phi)
        
        cpsi = cos(Bmag);			% cos(psi)
        spsi = sin(Bmag);			% sin(psi)
        
        if Bmag ~= 0
            Mx0 = M(lp-1,1);
            My0 = M(lp-1,2);
            Mz0 = M(lp-1,3);
            
            Mx1 = cphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
                + spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
                - sphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
                + cpsi.*(cphi.*My0-sphi.*Mx0));
            
            My1 = sphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
                + spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
                + cphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
                + cpsi.*(cphi.*My0-sphi.*Mx0));
            
            Mz1 = ct.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0)) ...
                - st.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
                + spsi.*(cphi.*My0-sphi.*Mx0));
        else
            Mx1 = M(lp-1,1);
            My1 = M(lp-1,2);
            Mz1 = M(lp-1,3);
        end
        
        % relaxation effects: "1" in Mz since Mo=1 by assumption
        M(lp,1) = Mx1 .* T2;
        M(lp,2) = My1 .* T2;
        M(lp,3) = Mz1 + (1 - Mz1) .* T1;
        
    end
    
end

