function WriteGenerationOutput_M300(fitness)

% USPEX Version 9.3.2
% Change: new output, hardness, order, etc

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resFolder = POP_STRUC.resFolder;

IND=POP_STRUC.ranking(1);
writeOUT_POSCAR(IND,1);
[nothing, nothing] = unix([' cat POSCAR >>      ' resFolder '/BESTgatheredPOSCARS']);
[nothing, nothing] = unix([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
update_USPEX_GENERATION(IND, fitness, 1);

% FIGURES
cd(resFolder);
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
extractStructures(20, ORG_STRUC.weight);
cd ..
WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);
safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
