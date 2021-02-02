function Start_POP_110()

global ORG_STRUC
global POP_STRUC

N_Step = length([ORG_STRUC.abinitioCode]);
POP_STRUC.DoneOrder = zeros(1, length(POP_STRUC.POPULATION));
for i = 1:length(POP_STRUC.POPULATION)
    if isempty(POP_STRUC.POPULATION(i).Step)
       POP_STRUC.POPULATION(i).Step = 1;
       POP_STRUC.POPULATION(i).Enthalpies = 100000*ones(1, N_Step);
       POP_STRUC.POPULATION(i).K_POINTS = ones(N_Step, 3);
    end
    POP_STRUC.POPULATION(i).Error = 0;
    POP_STRUC.POPULATION(i).Folder = 0;
    POP_STRUC.POPULATION(i).ToDo = 1;
    POP_STRUC.POPULATION(i).Done = 0;
    POP_STRUC.POPULATION(i).Number = 0;

    cattedCoors=[];
    typesAList = POP_STRUC.POPULATION(i).typesAList;
    numMols = POP_STRUC.POPULATION(i).numMols;
    for j = 1: sum(numMols)
        cattedCoors = cat(1,cattedCoors,POP_STRUC.POPULATION(i).MOLECULES(j).MOLCOORS);
    end
    POP_STRUC.POPULATION(i).COORDINATES = cattedCoors/(POP_STRUC.POPULATION(i).LATTICE);
    saveded = POP_STRUC.POPULATION(i).COORDINATES;
    newCoords = zeros(0,3);
    for m = 1: length(ORG_STRUC.atomType)
      s = find(typesAList == ORG_STRUC.atomType(m));
      newCoords = cat(1,newCoords, POP_STRUC.POPULATION(i).COORDINATES(s,:));
    end
    POP_STRUC.POPULATION(i).COORDINATES = newCoords;
    POP_STRUC.POPULATION(i).INIT_COORD = POP_STRUC.POPULATION(i).COORDINATES;
    POP_STRUC.POPULATION(i).INIT_LAT   = POP_STRUC.POPULATION(i).LATTICE;
    POP_STRUC.POPULATION(i).INIT_numIons=POP_STRUC.POPULATION(i).numIons;

end

