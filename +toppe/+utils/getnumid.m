function id = getNumId(sId)
% Convert three-letter string sId to numerical id for use with toppev3
% Use of this function is optional.

if length(sId) ~= 3 | ~ischar(sId)
   error('id must be a char array of length 3');
end

alph = 'abcdefghijklmnopqrstuvwxyz';
s = '';
for i = 1:length(sId)
	s = [s num2str(find(alph==sId(i)))];
end

id = str2num(s);

