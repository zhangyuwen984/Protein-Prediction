function getPOSCAR(ID, POSCARS_file)
% $Rev: 584 $
% $Author: maxim $
% $Date: 2014-08-30 10:27:47 +0400 (Sat, 30 Aug 2014) $

% Read all structures to the list:
%POSCARS_file = 'gatheredPOSCARS';
POSCAR_content = [];
fid   = fopen(POSCARS_file, 'rb');
tline = fgets(fid);
while ischar(tline)
    POSCAR_content = [POSCAR_content, tline];
    tline = fgets(fid);
end
fclose(fid);    % End of file reading

all_EA_occurrences = [strfind(POSCAR_content, 'EA'), length(POSCAR_content) + 1];

% Save found structures to temporary file:
f_poscar = fopen('POSCAR', 'wb');
EA_start = all_EA_occurrences(ID);
EA_next  = all_EA_occurrences(ID + 1);

structure = POSCAR_content(EA_start:EA_next - 1);
fprintf(f_poscar, structure);

fclose(f_poscar); % End of file writing

end