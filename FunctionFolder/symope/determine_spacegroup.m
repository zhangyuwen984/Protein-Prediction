function [resulted_sg_no, resulted_sg_name] = ...
    determine_spacegroup(lattice, coordinates, atom_type, ...
    num_atoms, nsym, structure_id, ...
    tolerance, home_path, method)
% The function determines symmetry using either Phonopy (if installed in the
% system) or Stokes' code (otherwise).
%
% Inputs:
%   lattice         : lattice;
%   coordinates     : coordinates;
%   atom_type       : types of atoms in the system;
%   num_atoms       : numbers of atoms;
%   nsym            : requested symmetry group;
%   structure_id    : a unique id to be added to the generated POSCAR name;
%   tolerance       : symmetry tolerance;
%   home_path       : a home path of USPEX calculation;
%   method          : optional parameter, can be 'phonopy' or 'stokes' (default).
% Outputs:
%   resulted_sg_no  : found space group number;
%   resulted_sg_name: found space group name.

if nargin < 9
    % Methods can be switched here:
    % method = 'Phonopy';
    method = 'Stokes';  % default method to determine space group symmetry
end


%% Parameters:
error_phonopy = 1;
if isequal(method, 'Phonopy')
    % Use Phonopy code (if it's found on a system) to determine the symmetry:
    [error_phonopy, nothing] = unix('which phonopy');
end

%% Determine real space group:
if error_phonopy == 0
    structures_dir = [home_path '/structures'];
    
    if ~isequal(exist(structures_dir, 'dir'), 7)  % 7 = directory
        [nothing, nothing] = unix(['mkdir ' structures_dir]);
    end
    
    structure_id = num2str(structure_id, '%04i');
    outfile = [structures_dir '/' structure_id '_SG_' num2str(nsym)];  % output POSCAR file
    
    
    %% Create POSCAR file - required for Phonopy and optional for Stokes code:
    atoms = '';
    for i=1:size(atom_type, 2)
        if i == 1
            atoms = [char(megaDoof(atom_type(i)))];
        else
            atoms = [atoms ' ' char(megaDoof(atom_type(i)))];
        end
    end
    atoms = {atoms};
    content = generate_poscar(lattice, coordinates, atoms, num_atoms, outfile, nsym);
    
    phonopy_command = ['ulimit -s unlimited; phonopy --symmetry --tolerance=' num2str(tolerance) ...
        ' -c ' outfile ' | grep space_group_type | cut -d: -f2- | awk ''{print $1,$2}'' '];
    [nothing, resulted_sg] = unix(phonopy_command);
    
    try
        parsed = sscanf(resulted_sg, '%s (%i)');
        
        resulted_sg_no   = parsed(end);
        resulted_sg_name = char(parsed(1:end-1)');
    catch
        resulted_sg_no   = 1;
        resulted_sg_name = 'P1';
    end
else
    % Use Stokes code to determine the symmetry:
    method = 'Stokes';
    resulted_sg_no = findsym_stokes(lattice, coordinates, num_atoms, atom_type, tolerance);
    if ~isempty(resulted_sg_no)
        resulted_sg_no = str2num(resulted_sg_no);
    else
        resulted_sg_no = 1;
    end
    resulted_sg_name = spaceGroups(resulted_sg_no);
end

%% Print report:
% disp(' ');
disp(['  Initial space group: ' spaceGroups(nsym) ' (' num2str(nsym)           ')']);
disp(['  Actual  space group: ' resulted_sg_name  ' (' num2str(resulted_sg_no) ...
    ') (determined with tolerance=' num2str(tolerance) ')']);
