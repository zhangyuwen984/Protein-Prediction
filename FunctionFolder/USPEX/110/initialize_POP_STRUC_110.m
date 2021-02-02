function initialize_POP_STRUC_110()

% USPEX Version 9.3.0
% Change: variable composition added
global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{},...
 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'convex_hull',{}, 'fitness', {});
POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {},'MOLECULES',{},...
  'S_order',{}, 'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {},'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
'struc_entr',{},'howCome',{},'JobID',{},'Folder',{}, 'numIons',{},'numMols',{},'MtypeLIST',{}, 'typesAList',{},'Number',{},'symg', {});
POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);
POP_STRUC(1).POPULATION(1).MOLECULES=struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{},'order',{});
POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'molecules',{},'fingerprint',{},'molfingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;
%create good initial population. Every individual fulfills hard constraints.
for i = 1: ORG_STRUC.initialPopSize
    [Molecules, lat, MtypeLIST, typesAList, numIons] = Random_Init_110(i, ORG_STRUC.numMols);
    POP_STRUC.POPULATION(i).LATTICE   = lat;
    POP_STRUC.POPULATION(i).MOLECULES = Molecules;
    POP_STRUC.POPULATION(i).MtypeLIST = MtypeLIST;
    POP_STRUC.POPULATION(i).typesAList= typesAList;
    POP_STRUC.POPULATION(i).numIons   = numIons;
    POP_STRUC.POPULATION(i).numMols   = ORG_STRUC.numMols;
    POP_STRUC.POPULATION(i).howCome = '  Random  ';
end

pick_Seeds();
Start_POP_110();
