function initialize_POP_STRUC_M200()
% 2D

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{},'SOFTMODEParents',{}, 'SOFTMUTATED',{},'resFolder', {},'generation',{},...
               'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'convex_hull',{}, 'fitness', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons',{},...
'struc_entr',{}, 'order',{},'FINGERPRINT',{},'K_POINTS',{},'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
'S_order',{},'howCome',{},'JobID',{},'Folder',{},'COORDINATES_2D',{},'LATTICE_2D',{},'Number',{}, 'symg',{},...
'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {});


POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

%create good initial population. Every individual fulfills hard constraints.
for i = 1: ORG_STRUC.initialPopSize
    if ~isempty(ORG_STRUC.numMols)
      [candidate_2D, lat_2D, candidate, lat] = Random_Init_M210(i, ORG_STRUC.numMols);
    else
      [candidate_2D, lat_2D, candidate, lat] = Random_Init_M200(i, ORG_STRUC.numIons);
    end
      POP_STRUC.POPULATION(i).LATTICE = lat;
      POP_STRUC.POPULATION(i).LATTICE_2D = lat_2D;
      POP_STRUC.POPULATION(i).COORDINATES = candidate;
      POP_STRUC.POPULATION(i).COORDINATES_2D = candidate_2D;
      POP_STRUC.POPULATION(i).howCome = '  Random  ';
      POP_STRUC.POPULATION(i).numIons = ORG_STRUC.numIons;
end

if ORG_STRUC.spin == 1
    for i = 1: ORG_STRUC.initialPopSize
        POP_STRUC.POPULATION(i) = individual_Spin_Init(  POP_STRUC.POPULATION(i) );
    end
end

    pick_Seeds();
Start_POP_M200();
