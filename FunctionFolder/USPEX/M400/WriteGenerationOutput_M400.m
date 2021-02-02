function WriteGenerationOutput_M400(fitness)
% $Rev: 1120 $
% $Author: mrakitin $
% $Date: 2015-08-20 21:51:42 +0400 (Thu, 20 Aug 2015) $

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resFolder = POP_STRUC.resFolder;

IND=POP_STRUC.ranking(1);

update_USPEX_GENERATION(IND, fitness, 1);

cd(resFolder);

getPDB(POP_STRUC.POPULATION(IND).Number, 'gatheredPDB');
[nothing, nothing] = unix(['cat PDB >> BESTgatheredPDB']);

pdb_folder = 'BEST_PDB';
if ~isequal(exist(pdb_folder, 'dir'), 7) % 7 = directory.
    mkdir(pdb_folder);
end
[nothing, nothing] = unix(['cp -f PDB ' pdb_folder '/EA' sprintf('%04d', POP_STRUC.POPULATION(IND).Number) '.pdb']);
[nothing, nothing] = unix(['rm -f PDB']);

getPOSCAR(POP_STRUC.POPULATION(IND).Number, 'gatheredPOSCARS');
[nothing, nothing] = unix(['cat POSCAR >> BESTgatheredPOSCARS']);
[nothing, nothing] = unix(['rm -f POSCAR']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RMSD OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fit_loop = 1:length(POP_STRUC.POPULATION)
    if ~isempty(POP_STRUC.POPULATION(fit_loop).RMSD)	
        fitness_rmsd(fit_loop) = POP_STRUC.POPULATION(fit_loop).RMSD(end);
    else
        fitness_rmsd(fit_loop) = 10000
    end
end
[nothing, ranking_rmsd] = sort(fitness_rmsd);
IND=ranking_rmsd(1);
getPDB(POP_STRUC.POPULATION(IND).Number, 'gatheredPDB');

pdb_folder = 'BEST_RMSD';
if ~isequal(exist(pdb_folder, 'dir'), 7) % 7 = directory.
    mkdir(pdb_folder);
end
[nothing, nothing] = unix(['cp -f PDB ' pdb_folder '/EA' sprintf('%04d', POP_STRUC.POPULATION(IND).Number) '.pdb']);
[nothing, nothing] = unix(['rm -f PDB']);


% FIGURES
makeFigures_M400(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
cd ..

WriteIndividual(resFolder, 'kcal/mol');
WriteBest(resFolder, 'kcal/mol');
WriteGeneration(resFolder);

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
