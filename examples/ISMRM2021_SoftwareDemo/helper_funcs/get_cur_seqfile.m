function [seqfile_str] = get_cur_seqfile()
% get the "current" .seq file, which is defined as the most recently
% written .seq file in the directory, and return it as a string

S = dir('*.seq');
S = S(~[S.isdir]);
[~,idx] = sort([S.datenum]);
S = S(idx);
S_last = S(end);
seqfile_str = S_last.name;

end

