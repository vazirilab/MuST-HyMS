%%   Get value from ScanImage meta data.
% 
% See github wiki for documentation. 
%
% Dependencies:
% -
% 
% Input
% metadata: Meta data from ScanImage TIF file.
% str: String for meta data variable.
% num_length: Length of the number.
%
% Output
% out: Value for the meta data variable.
%
%
% --SW, last modified: 12/14/2018.


function [out] = GetMetaDataNumber(metadata,str,num_length)

ind = strfind(metadata,str);
start_ind = ind + length(str) + 3;
end_ind = start_ind + num_length;
str_out = metadata(start_ind:end_ind);
str_number = regexp(str_out,'(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?','match');
out = str2double(str_number);

end