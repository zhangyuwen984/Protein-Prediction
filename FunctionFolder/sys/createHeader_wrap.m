function formatted_rows = createHeader_wrap(text, position)

if nargin < 2
    position = 'center';
end

width    = 80;

%---> Format text:
line_len       = width;
formatted_rows = {};
delimiter      = '|';
delimiter2     = '-';
separator      = ['*' strrep(blanks(line_len-2), ' ', delimiter2) '*'];

% Beginning separator:
formatted_rows{end+1} = separator;


blanks_before = 2;
for i=1:size(text,2)
    if strcmp(position, 'center') > 0
        blanks_before = round(line_len/2 - size(text{i},2)/2 - 1);
    end
    blanks_after = line_len - (1 + blanks_before + size(text{i},2)) - 1;
    a = [delimiter blanks(blanks_before) text{i} blanks(blanks_after) delimiter];
    
    formatted_rows{end+1} = a;
end

% Ending separator:
formatted_rows{end+1} = separator;

end
