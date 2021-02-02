function PROTEINS_STRUC = backbone_int()
% $Rev: 1120 $
% $Author: mrakitin $
% $Date: 2015-08-20 21:51:42 +0400 (Thu, 20 Aug 2015) $

% This function gets information about protein backbone from input.make file.

global ORG_STRUC

PROTEINS_STRUC = struct(...
                        'residue_rows'   , '', ...
                        'residue_numbers', '', ...
                        'backbone_atoms' , '', ...
                        'backbone_str'   , '' ...
                       );

% Temporary directory for backbone calculation:
backbone_dir = 'backboneCalc';

% Clean before execution:
if isequal(exist(backbone_dir, 'dir'), 7) % 7 = directory
    rmdir(backbone_dir, 's');
end
mkdir(backbone_dir);

% Copy all files from Specific/ directory to temp directory to execute
% protein.x and get input.int (internal coordinates with marker angles 
% 0.001, 0.002, etc.):
[nothing, nothing] = unix(['cp ./'  ORG_STRUC.specificFolder '/* ' backbone_dir]);
cd(backbone_dir);
copyfile('input.make0_1', 'input.make', 'f');
copyfile('input.key_1'  , 'input.key' , 'f');


% Read input.make input file and read lines and numbers of residues:
make_content = {};
fid = fopen('input.make', 'rb');
tline = fgets(fid);
i = 0;
while ischar(tline)
    i = i + 1;
    make_content{i} = tline;
    tline = fgets(fid);
end
fclose(fid);    % End of file reading

residue_start = -1;
residue_end   = -1;

for i=1:size(make_content,2)
    current_str = make_content{i};
    search = strfind(current_str, '0.001');
    if ~isempty(search)
        residue_start = i;
        break
    end
end

counter = 0;
rows    = [];
values  = [];
for i=residue_start:size(make_content,2)
    counter = counter + 1;

    current_str = make_content{i};
    
    %splitted    = strsplit(current_str);    
    % strsplit doesn't exist in older versions of Matlab and appeared
    % recently in Matlab R2013a+ version:
    % http://stackoverflow.com/questions/18673969/strsplit-undefined-function-for-input-type-char?answertab=active#tab-top
    % Since textscan will be implemented as a workaround:
    
    splitted    = textscan(current_str, '%s %f %f %f');
    value_x1000 = splitted{2}*1000;
    if ~isequal(counter, value_x1000)
        break
    end
    residue_end = i;
    rows   = [rows, i];
    values = [values, value_x1000];
end


% Read numbers of C atoms in input.int:
[nothing, nothing] = unix('protein < input.make > /dev/null');

int_content = {};
fid = fopen('input.int', 'rb');
tline = fgets(fid);
i = 0;
while ischar(tline)
    i = i + 1;
    int_content{i} = tline;
    tline = fgets(fid);
end
fclose(fid);    % End of file reading

backbone_str = {};
backbone_num = [];
counter = 0;
for i=1:size(int_content,2)
    current_str = int_content{i};
    for j=1:size(values,2)
        search = strfind(current_str, [' ' sprintf('%0.4f',(values(j)/1000)) ' ']);
        if ~isempty(search)
            C_num = textscan(current_str, '%s');

            backbone_num = [backbone_num, str2num(C_num{1}{1})];

            counter = counter + 1;
            backbone_str{counter} = current_str;
        end
    end
end

backbone_num = [3, backbone_num(2:end)]

% Create return struct:
PROTEINS_STRUC.residue_rows    = rows;
PROTEINS_STRUC.residue_numbers = values;
PROTEINS_STRUC.backbone_atoms  = backbone_num
PROTEINS_STRUC.backbone_str    = backbone_str;


cd('..');

% Clean after execution:
%[nothing, nothing] = unix(['rm -rf ' backbone_dir '/']);
if isequal(exist(backbone_dir, 'dir'), 7) % 7 = directory
    rmdir(backbone_dir, 's');
end

end
