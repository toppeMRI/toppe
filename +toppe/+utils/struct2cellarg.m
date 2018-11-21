function c = struct2cellarg(s)
% Convert struct into paired cell array {'field1', value1, 'field2', value2, ...}.
%
% function ca = struct2cellarg(s)
%
% Basically does the opposite of built-in struct() function, which converts
% paired cell array {'field1', 'value1', 'field2', 'value2', ...} into a struct.
%
% Useful for passing paired varargin-type arguments to functions, 
% by passing as {:} (comma-separated list)
%
% Input:
%   s     struct 
%           s.field1 = value1;
%           s.field2 = value2; 
%           etc
% Output:
%   c     paired cell array: {'field1', value1, 'field2', value2, ...}
%
% Example usage:
%  >> opts.rf = rf;
%  >> opts.gz = gz;
%  >> opts.ofname = 'mytoppemodule.mod';
%  >> opts = struct2cellarg(opts);
%  >> writemod(opts{:});
%
% $Id: struct2cellarg.m,v 1.2 2018/10/24 10:44:50 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/struct2cellarg.m,v $

fields = fieldnames(s);
vals = struct2cell(s);
c = [];
for ii = 1:length(fields)
	c = [c fields(ii) vals(ii)];
end
