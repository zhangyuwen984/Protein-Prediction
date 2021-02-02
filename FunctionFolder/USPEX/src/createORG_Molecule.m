function createORG_Molecule(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%[nothing, checkMolecules] =unix (['./getStuff ' inputFile ' checkMolecules 1']);
checkMolecules = python_uspex(getPy, ['-f ' inputFile ' -b checkMolecules -c 1']);
if ~isempty(checkMolecules)
    ORG_STRUC.checkMolecules = str2num(checkMolecules);
end
 
  
read_molecules();
