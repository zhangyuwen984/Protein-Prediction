function Start_POP_201()

global ORG_STRUC
global POP_STRUC

N_Step = length([ORG_STRUC.abinitioCode]);
POP_STRUC.DoneOrder = zeros(1, length(POP_STRUC.POPULATION));
for i = 1:length(POP_STRUC.POPULATION)
    POP_STRUC.POPULATION(i).INIT_COORD = POP_STRUC.POPULATION(i).COORDINATES;
    POP_STRUC.POPULATION(i).INIT_LAT   = POP_STRUC.POPULATION(i).LATTICE;
    POP_STRUC.POPULATION(i).INIT_numIons=POP_STRUC.POPULATION(i).numIons;
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
end

ReRank();
