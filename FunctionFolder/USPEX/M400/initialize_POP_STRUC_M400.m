function initialize_POP_STRUC_M400()
% $Rev: 690 $
% $Author: maxim $
% $Date: 2014-10-28 07:25:19 +0400 (Tue, 28 Oct 2014) $

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{}, ...
                  'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'fitness', {}, 'convex_hull',{});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'numIons', {}, 'INIT_numIons', {},...
'struc_entr',{}, 'order',{},'dielectric_tensor',{}, 'gap',{}, 'hardness',{}, 'mag_moment',{}, 'FINGERPRINT',{}, 'symg',{}, 'K_POINTS',{},...
'S_order',{}, 'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},'howCome',{}, 'template_random', {}, 'JobID',{},'Folder',{}, 'Number',{});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount  = 0;
POP_STRUC.bad_rank   = 0;

% Get number of backbone atoms:
PROTEINS_STRUC = backbone_int();
POP_STRUC.backbone_atoms = PROTEINS_STRUC.backbone_atoms;

% Create a good initial population:
for i = 1: ORG_STRUC.initialPopSize 

    angles_num         = size(POP_STRUC.backbone_atoms, 2);

    %POP_STRUC.POPULATION(i).LATTICE     = '';
    %POP_STRUC.POPULATION(i).COORDINATES = '';

    [POP_STRUC.POPULATION(i).ANGLES, POP_STRUC.POPULATION(i).template_random]  = Random_Init_M400(angles_num);
    POP_STRUC.POPULATION(i).howCome = 'Random';
    POP_STRUC.POPULATION(i).numIons = angles_num;
end

%%%%%%%%%%%% END defining the first generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% SEEDING %%%%%%%%%%%%%%%%%
pick_Seeds();

if ORG_STRUC.doFing
    pickAntiSeeds();
end

Start_POP_M400();
