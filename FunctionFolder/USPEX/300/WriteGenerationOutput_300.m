function WriteGenerationOutput_300(fitness)

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;

IND=POP_STRUC.ranking(1);
lattice = POP_STRUC.POPULATION(IND).LATTICE;
coor    = POP_STRUC.POPULATION(IND).COORDINATES;
numIons = POP_STRUC.POPULATION(IND).numIons;
lattice = POP_STRUC.POPULATION(IND).LATTICE;
symg    = POP_STRUC.POPULATION(IND).symg;
order   = POP_STRUC.POPULATION(IND).order;
count   = POP_STRUC.POPULATION(IND).Number;

Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
[nothing, nothing] = unix([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
[nothing, nothing] = unix([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);

update_USPEX_GENERATION(IND, fitness, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make Sure bodyCount is correct

% FIGURES
cd (resFolder);
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
extractStructures(ORG_STRUC.populationSize, ORG_STRUC.weight);
if ORG_STRUC.PhaseDiagram
    phase_diagram;
end

V = [];
for i=1:length(USPEX_STRUC.POPULATION)
    V = [V det(USPEX_STRUC.POPULATION(i).LATTICE)*sum(ORG_STRUC.numIons)/sum(USPEX_STRUC.POPULATION(i).numIons) ];
end
ORG_STRUC.latVolume = mean(V);

cd ..

fpath = [ resFolder '/OUTPUT.txt'];
fp = fopen(fpath, 'a+');
fprintf(fp, [alignLine( sprintf('Approximate volume for each species: %.4f A^3', mean(V)) ) '\n'] );
fclose(fp);

WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);
WriteProperties(ORG_STRUC.optType, resFolder);
WriteCompStatistic(resFolder)

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
