function blk = readblock(fname)

C = toppe.constants;

fid = fopen(fname, 'r', 'ieee-be');

% rf
type = fread(fid, 1, 'int16');
if type ~= C.NULL
    blk.rf.delay = fread(fid, 1, 'int16')/1e6;  % sec
    blk.rf.amp = fread(fid, 1, 'int16')/C.RFSCALE;     % Gauss
    n = fread(fid, 1, 'int16');
    blk.rf.rho = fread(fid, n, 'int16');
    blk.rf.theta = fread(fid, n, 'int16');
    %rho = fread(fid, n, 'int16')/C.MAXIAMP*blk.rf.amp;    % Gauss
    %theta = fread(fid, n, 'int16')/C.MAXIAMP*pi;   % radians
    %blk.rf.signal = rho.*exp(1i*theta);
else
    blk.rf = [];
end

% gradients 
blk.gx = sub_readgrad(fid);
blk.gy = sub_readgrad(fid);
blk.gz = sub_readgrad(fid);

% adc
type = fread(fid, 1, 'int16');
if type == C.ADC
    blk.adc.numSamples = fread(fid, 1, 'int16');
    blk.adc.dwell = fread(fid, 1, 'int16')/1e6;
    blk.adc.delay = fread(fid, 1, 'int16')/1e6;
else
    blk.adc = [];
end

fclose(fid);

return


function g = sub_readgrad(fid)

C = toppe.constants;

type = fread(fid, 1, 'int16');

if type == C.NULL
    g = [];
    return;
end

delayus = fread(fid, 1, 'int16');
g.delay = delayus/1e6;  % sec
g.amplitude = fread(fid, 1, 'int16')/C.GSCALE;   % Gauss/cm

if type == C.TRAP
    g.type = 'trap';
    g.riseTime = fread(fid, 1, 'int16')/1e6;
    g.flatTime = fread(fid, 1, 'int16')/1e6;
    g.fallTime = fread(fid, 1, 'int16')/1e6;
else
    % TODO
end

return
