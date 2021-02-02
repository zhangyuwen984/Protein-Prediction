function xyz_backbone = backbone_xyz(backbone_atoms, xyzfile)
% $Rev: 584 $
% $Author: maxim $
% $Date: 2014-08-30 10:27:47 +0400 (Sat, 30 Aug 2014) $

% Read coordinates from input.xyz:
xyz_content = {};
fid = fopen(xyzfile, 'rb');
tline = fgets(fid);
i = 0;
while ischar(tline)
    i = i + 1;
    xyz_content{i} = tline;
    tline = fgets(fid);
end
fclose(fid);    % End of file reading


xyz_str      = {};
xyz_backbone = [];
counter = 0;
for i=1:size(xyz_content,2)
    current_str = xyz_content{i};
    for j=1:size(backbone_atoms,2)
        splitted = textscan(current_str, '%s');
        if strcmp(splitted{1}{1}, num2str(backbone_atoms(j))) == 1 & strcmp(splitted{1}{2}, 'C') == 1
            counter = counter + 1;
            xyz_str{counter} = current_str;            
            xyz_backbone = [xyz_backbone; str2num(splitted{1}{3}) str2num(splitted{1}{4}) str2num(splitted{1}{5})];
        end
    end
end

end
