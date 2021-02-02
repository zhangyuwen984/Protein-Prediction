function lineout = alignLine(linein, blank_surround)

if nargin < 2
    blank_surround = 1;
end

%---> Format text:

width          = 80;
line_len       = width;
lineout        = linein;
delimiter      = '-';
blanks_before  = round(line_len/2 - size(linein,2)/2);
blanks_after   = line_len - (blanks_before + size(linein,2));
if blank_surround == 1
    space = ' ';
else
    space = delimiter;    
end

lineout        = [strrep(blanks(blanks_before-1), ' ', delimiter) space linein space strrep(blanks(blanks_after-1), ' ', delimiter)];

end
