function writecoresfile(cores)

fid = fopen('cores.txt', 'wt');
fprintf(fid, 'Total number of cores\n');
fprintf(fid, '%d\n', length(cores));
fprintf(fid, 'nmodules modIds... \n');

% Recall that module id's are defined by their row number in modules.txt
for icore = 1:length(cores)
    nmod = length(cores{icore});
    fprintf(fid, '%d\t', nmod);
    for imod = 1:nmod
        modid = cores{icore}(imod);
        fprintf(fid, '%d', modid);
        if imod < nmod
            fprintf(fid, '\t');
        else
            fprintf(fid, '\n');
        end
    end
end
