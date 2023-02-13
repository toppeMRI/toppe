function blk = readblock(fname)

C = toppe.constants;

fid = fopen(fname, 'r', 'ieee-be');

% rf
type = fread(fid, 1, 'int16');
if type ~= C.NULL
    delayus = fread(fid, 1, 'int16');
    blk.rf.delay = delayus/1e6;  % sec
    amp = fread(fid, 1, 'int16')/C.rfscale;
    n = fread(fid, 1, 'int16');
    rho = fread(fid, n, 'int16')/C.max_pg_iamp*amp;  % Gauss
    theta = fread(fid, n, 'int16');
    blk.rf.signal = rho.*exp(1i*theta/C.max_pg_iamp);
else
    blk.rf = [];
end

% gradients (remember that order is important)
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
amp = fread(fid, 1, 'int16')/C.gscale;

if type == C.TRAP
    g.type = 'trap';
    g.riseTime = fread(fid, 1, 'int16')/1e6;
    g.flatTime = fread(fid, 1, 'int16')/1e6;
    g.fallTime = fread(fid, 1, 'int16')/1e6;
else
    % TODO
end

return
