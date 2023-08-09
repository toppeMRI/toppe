function [pfile_str] = get_cur_pfile()
% get the "current" P-file, which is defined as the most recently
% written .7 file in the directory:
%
% /mnt/storage/mhaskell/data/yyyy/mm/yyyy_mm_dd
%
% andreturn it as a string

%% Set the path to look for the P-file data
m = datestr(now,'mm'); y = datestr(now,'yyyy'); 
path = ['/mnt/storage/mhaskell/data/',y,'/',m,'/',datestr(now,'yyyy_mm_dd'),'/'];


%% Go to the data directory to get data info, then return to current dir
cur_dir = pwd;
cd(path)
try
    S = dir('*.7');
catch
    disp('error getting directory info')
end
cd(cur_dir)


%% Find the most recent P-file written and return its full filepath
S = S(~[S.isdir]);
[~,idx] = sort([S.datenum]);
S = S(idx);
S_last = S(end);
pfile_str = [path, S_last.name];


end

