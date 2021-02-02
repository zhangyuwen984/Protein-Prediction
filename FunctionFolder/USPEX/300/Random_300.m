function Random_300(Ind_No)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global OFF_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSplit    = length( ORG_STRUC.firstGeneSplit );
numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );

k=1;
while ~CompositionCheck( numBlocks )
    numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
    k = k+1;
    if k > 1000
        disp('Error in Anti-compositions file, please check !!');
        quit
    end
end

numIons = ORG_STRUC.numIons*numBlocks;
numMols = ORG_STRUC.numMols*numBlocks;
if ~isempty(ORG_STRUC.numMols)
    [candidate, lat] = Random_Init_310(Ind_No, numMols);
else
    [candidate, lat] = Random_Init_300(Ind_No, numIons);
end

OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
OFF_STRUC.POPULATION(Ind_No).numIons = numIons;