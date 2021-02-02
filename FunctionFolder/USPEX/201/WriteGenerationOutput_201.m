function WriteGenerationOutput_201(fitness)

% USPEX Version 9.3.2

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;
convex_hull = POP_STRUC.convex_hull;

IND=[];
if ~isempty(convex_hull)
    fpath1 = [ resFolder '/convex_hull'];
    fp1 = fopen(fpath1, 'a+');
    fprintf(fp1, '---- generation%3d ----\n', POP_STRUC.generation);
    item = 1;
    for i = 1:size(convex_hull,1)
        if convex_hull(i,3) > 0 % stable composition
            IND(item) = round(convex_hull(i,3));
            fprintf(fp1,'%8.4f  %12.4f\n', convex_hull(i,1:2));
            item = item + 1;
        end
    end
    fclose(fp1);
else
    IND=POP_STRUC.ranking(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make Sure bodyCount is correct
for i=1:length(IND)
    lattice = POP_STRUC.POPULATION(IND(i)).LATTICE;
    coor    = POP_STRUC.POPULATION(IND(i)).COORDINATES;
    numIons = POP_STRUC.POPULATION(IND(i)).numIons;
    lattice = POP_STRUC.POPULATION(IND(i)).LATTICE;
    symg    = POP_STRUC.POPULATION(IND(i)).symg;
    %    order   = POP_STRUC.POPULATION(IND(i)).order;
    
    count1 = 0;
    order = zeros(sum(numIons),1);
    Surface_order   = POP_STRUC.POPULATION(IND(i)).Surface_order;
    for j=1:size(coor,1)
        if POP_STRUC.POPULATION(IND(i)).chanAList(j)==1
            count1=count1+1;
            order(j) = Surface_order(count1);
        end
    end
    
    count   = POP_STRUC.POPULATION(IND(i)).Number;
    Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
    Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
    
    [nothing, nothing] = unix([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
    [nothing, nothing] = unix([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
end
update_USPEX_GENERATION(IND, fitness, 2);
% FIGURES

cd(resFolder);
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);

if sum(ORG_STRUC.bulk_ntyp >0) == 1   %substrate is elemental
    if length(ORG_STRUC.numIons) == 2 % only one type of surface atom(eg: C on C(111))
        extendedConvexHull_201(convex_hull, ORG_STRUC.E_A, ORG_STRUC.bulk_stoi, 1.0, 1);
    else
        extractStructures(20, ORG_STRUC.weight);
    end
else
    extendedConvexHull_201(convex_hull, ORG_STRUC.E_AB,ORG_STRUC.bulk_stoi, 1.0, 2);
end
%     extractDifferentStructures_var(0.30, 0.04, sum(ORG_STRUC.numIons));
cd ..
WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);

if exist('writeMagmoment')
    writeMagmoment(resFolder);
end

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);


