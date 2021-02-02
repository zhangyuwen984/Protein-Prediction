function formatted_rows = createHeader(pic_file, description, funcFold)
% The function is to create headers for USPEX, META, etc. in OUTPUT.txt

pic_file = [funcFold '/' pic_file];
versfile = [funcFold '/VERSION'];
version  = '';
fid      = fopen(versfile, 'rb');
version  = deblank(fgets(fid));
fclose(fid);

width    = 80;
text     = {version, ...
            '', ...
            description, ...
            'more info at http://uspex.stonybrook.edu', ...
            '' ...
            };

fid     = fopen(pic_file, 'rb');
tline   = fgets(fid);
rows    = {};
max_len = 0;
while ischar(tline)
    tline       = deblank(tline);
    rows{end+1} = tline;
    new_max_len = size(tline,2);
    if new_max_len > max_len
        max_len = new_max_len;
    end
    tline = fgets(fid);
end
fclose(fid);

%---> Format picture:
line_len       = width;
formatted_rows = {};
delimiter      = '|';
delimiter2     = '-';
separator      = ['*' strrep(blanks(line_len-2), ' ', delimiter2) '*'];

% Beginning separator:
%disp(separator);
formatted_rows{end+1} = separator;

blanks_before = round(line_len/2 - max_len/2 - 1);
blanks_after  = 0;
for i=1:size(rows,2)
    blanks_after = line_len - (1 + blanks_before + size(rows{i},2)) - 1;
    a = [delimiter blanks(blanks_before) rows{i} blanks(blanks_after) delimiter];

    %disp(a);
    formatted_rows{end+1} = a;
end


%---> Format text:
for i=1:size(text,2)
    blanks_before = round(line_len/2 - size(text{i},2)/2 - 1);
    blanks_after = line_len - (1 + blanks_before + size(text{i},2)) - 1;
    a = [delimiter blanks(blanks_before) text{i} blanks(blanks_after) delimiter];

    %disp(a);
    formatted_rows{end+1} = a;
end


% Ending separator:
%disp(separator);
formatted_rows{end+1} = separator;

end
