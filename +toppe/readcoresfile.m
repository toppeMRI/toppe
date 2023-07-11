function blockGroups = readcoresfile(fn)

fid = fopen(fn, 'r');

s = fgets(fid);  % skip line
nGroups = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line

% Recall that module id's are defined by their row number in modules.txt
for i = 1:nGroups
    nBlocksInGroup = fscanf(fid, '%d\n', 1);
    clear blockIDs;
    for j = 1:nBlocksInGroup
        if j < nBlocksInGroup
            blockIDs(j) = fscanf(fid, '%d\t', 1);
        else
            blockIDs(j) = fscanf(fid, '%d\n', 1);
        end
    end
    blockGroups{i} = blockIDs;
end

fclose(fid);
