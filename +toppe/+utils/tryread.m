function out = tryread(fn, arg)
% Wrap function call in try/catch block
try
	out = fn(arg);
catch ME
	error('Failed to read %s\n', arg);
end
