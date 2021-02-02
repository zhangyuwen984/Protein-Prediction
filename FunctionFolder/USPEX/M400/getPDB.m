function getPDB(ID, PDB_file)
% $Rev: 637 $
% $Author: maxim $
% $Date: 2014-10-03 11:41:44 +0400 (Fri, 03 Oct 2014) $

% Read all structures to the list:
%PDB_file = 'gatheredPDB';
PDB_content = [];
fid   = fopen(PDB_file, 'rb');
tline = fgets(fid);
while ischar(tline)
    PDB_content = [PDB_content, tline];
    tline = fgets(fid);
end
fclose(fid);    % End of file reading

all_EA_occurrences = [strfind(PDB_content, 'HEADER    EA'), length(PDB_content) + 1];
current_EA         = [strfind(PDB_content, ['HEADER    EA' num2str(ID) ' ']), length(PDB_content) + 1];
struct_index       = find(all_EA_occurrences == current_EA(1));

% Save found structures to temporary file:
f_pdb = fopen('PDB', 'wb');
EA_start = all_EA_occurrences(struct_index);
EA_next  = all_EA_occurrences(struct_index + 1);

structure = PDB_content(EA_start:EA_next - 1);
fprintf(f_pdb, structure);

fclose(f_pdb); % End of file writing

end
