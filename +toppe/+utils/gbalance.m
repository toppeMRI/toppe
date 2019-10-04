function gbal = gbalance(g,smax,targetlength)

% G is input gradient waveform to balance
% smax is maximum slew rate
% targetlength produces an optimized waveform of specified length
% Typically you can balance Gx and Gy and one ends up being longer,
% then rerun the optimization on the shorter one so it can use the extra
% time

% Returns gradient waveform that balances the input
gmax = 5;
%smax = 20; %G/cm/ms
dt = 4e-3; %ms


%% Create initial solution which slews to 0 then makes a blip
rampsamps = ceil(abs(g(end)/(smax*dt)));
postramp = linspace(g(end),0,rampsamps);

m0 = sum([g; postramp'])*dt/1000;
if m0 < 0    
    blip = -toppe.utils.trapwave2(-m0,gmax,smax, 4e-3);
else
    blip = toppe.utils.trapwave2(m0,gmax,smax, 4e-3);
end
blip = blip(find(blip,1,'first'):end); % Remove leading zeros
x0 = [postramp'; -blip'];
totallength = length(g) + length(x0);
npad = 4-mod(totallength,4);
x0 = [x0; zeros(npad,1)];
x0length = length(x0);

if targetlength > 0 % Stretch x0 to match extra time
    x0 = interp1(x0,linspace(1,x0length,targetlength)).'/(targetlength/x0length);
end 

N = length(x0);

Aeq = [ones(1,N); zeros(1,N-1) 1];
Beq = [-trapz(g); 0];

% fmincon for flexibility even though this is a quadratic problem
options = optimoptions('fmincon','MaxFunctionEvaluations',inf,'Display','notify','Algorithm','sqp');
gbal = fmincon(@gbalcost,x0,[],[],Aeq,Beq,-gmax*ones(N,1),gmax*ones(N,1),@syscon,options);

% Check if ending is close enough to zero to round down
% Gmax = 5 / 32766 DAQ steps = 1.525e-4
if gbal(end) < 1.5e-4
    gbal(end) = 0;
else
    error('Failed to make ending value zero');
end

return

    function cost = gbalcost(x)
        % Total variation of slew rate, affects slew acc
        % Also account for the last sample of g
        tv = sum(abs(diff([g(end-1:end); x],2).^2));
        slew = sum(abs(diff([g(end-1:end); x]).^2));
        
        % Regularization parameters, affects smoothness
        lam = 1; % Acc
        lam2 = 0.02; % Vel
       
        
        cost = lam*tv + lam2*slew;
    end

    function [c,ceq] = syscon(x)
        % Hardware constraints
        
        slew = abs(diff([g(end); x]))/dt;
        
        c = slew - smax; % solver contraint is c <= 0
        ceq = zeros(size(x));
    end
end
