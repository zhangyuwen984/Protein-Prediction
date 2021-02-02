function initialize_POP_STRUC_310()

% USPEX Version 9.3.0
% Change: variable composition added
global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{},...
    'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'convex_hull',{}, 'fitness', {});

POP_STRUC(1).POPULATION = struct('COORDINATES',{}, 'INIT_COORD',{}, 'LATTICE',{}, 'INIT_LAT',{}, 'numIons',{}, 'INIT_numIons',{}, ...
    'struc_entr',{},'MOLECULES',{},'order',{},'FINGERPRINT',{},'K_POINTS',{},'Step',{},'Enthalpies', {},'Error',{},'Done',{},'ToDo',{}, ...
    'S_order',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{},'numMols',{},'MtypeLIST',{},'typesAList',{},'softmode',{},'Number',{},'symg',{});
POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).POPULATION(1).MOLECULES=struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{},'order',{},'Operation',{},'P',{},'PB',{});
POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'molecules',{},'fingerprint',{},'eignFre',{},...
    'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;


% create the composition file for 310
createCompostion_SingleBlock();


%create good initial population. Every individual fulfills hard constraints.
successPop = 0;
goodPop    = 0;
nsym       = 0;

nSplit = length( ORG_STRUC.firstGeneSplit );
if nSplit == 0
    error('ERROR:  The minAt and maxAt is too small, please check your input file ...');
end

numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
numMols = ORG_STRUC.numMols*numBlocks;

while successPop < ORG_STRUC.initialPopSize

    [goodPop, POP, nsym] = Random_Init_310(goodPop, numMols, nsym);
    if goodPop > successPop
        POP_STRUC.POPULATION(goodPop).MOLECULES  = POP.MOLECULES;
        POP_STRUC.POPULATION(goodPop).numMols    = POP.numMols;
        POP_STRUC.POPULATION(goodPop).MtypeLIST  = POP.MtypeLIST;
        POP_STRUC.POPULATION(goodPop).typesAList = POP.typesAList;
        POP_STRUC.POPULATION(goodPop).numIons    = POP.numIons;
        POP_STRUC.POPULATION(goodPop).LATTICE    = POP.LATTICE;
        POP_STRUC.POPULATION(goodPop).howCome    = '  Random  ';
        successPop = successPop +1;

        % To balance the compositions here !
        numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
        numMols = ORG_STRUC.numMols*numBlocks;
    end
end


pick_Seeds();
Start_POP_310();
%======================= END defining the first generation  =======================
%==================================================================================

%% This inline function is created to

function createCompostion_SingleBlock()

global ORG_STRUC


if isempty(ORG_STRUC.minAt) | isempty(ORG_STRUC.maxAt)
    ORG_STRUC.minAt = sum(ORG_STRUC.numIons);
    ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end

N_T = size(ORG_STRUC.numIons,1);
splitting = zeros(1,N_T);
findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
IPS = ORG_STRUC.initialPopSize;

fp = fopen('Seeds/compositions', 'w');
if exist('Seeds/Anti-compositions')
    [nothing, nothing] = unix('mv Seeds/Anti-compositions Seeds/Anti-compositions-back');
end
for i=1:size(ORG_STRUC.firstGeneSplit,1)
    for j=1:size(ORG_STRUC.firstGeneSplit,2)
        fprintf(fp, '%4d', ORG_STRUC.firstGeneSplit(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);
