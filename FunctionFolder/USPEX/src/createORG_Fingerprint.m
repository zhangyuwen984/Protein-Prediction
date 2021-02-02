function createORG_Fingerprint(inputFile)

global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fingerprints:
%[nothing, RmaxFing] = unix (['./getStuff ' inputFile ' RmaxFing 1']);
RmaxFing =  python_uspex(getPy, ['-f ' inputFile ' -b RmaxFing -c 1']);
if ~isempty(RmaxFing)
    ORG_STRUC.RmaxFing = str2num(RmaxFing);
end
if ORG_STRUC.RmaxFing > 0
    ORG_STRUC.doFing = 1;
else
    ORG_STRUC.doFing = 0;
end

%[nothing, deltaFing] = unix (['./getStuff ' inputFile ' deltaFing 1']);
deltaFing = python_uspex(getPy, ['-f ' inputFile ' -b deltaFing -c 1']);
if ~isempty(deltaFing)
    ORG_STRUC.deltaFing = str2num(deltaFing);
end

%[nothing, sigmaFing] = unix (['./getStuff ' inputFile ' sigmaFing 1']);
sigmaFing = python_uspex(getPy, ['-f ' inputFile ' -b sigmaFing -c 1']);
if ~isempty(sigmaFing)
    ORG_STRUC.sigmaFing = str2num(sigmaFing);
end

%[nothing, toleranceFing] = unix (['./getStuff ' inputFile ' toleranceFing 1']);
toleranceFing = python_uspex(getPy, ['-f ' inputFile ' -b toleranceFing -c 1']);
if ~isempty(toleranceFing)
    ORG_STRUC.toleranceFing = str2num(toleranceFing);
elseif ORG_STRUC.molecule
    ORG_STRUC.toleranceFing = 0.05;
end
%[nothing, toleranceBestHM] = unix (['./getStuff ' inputFile ' toleranceBestHM 1']);
toleranceBestHM = python_uspex(getPy, ['-f ' inputFile ' -b toleranceBestHM -c 1']);
if ~isempty(toleranceBestHM)
    ORG_STRUC.toleranceBestHM = str2num(toleranceBestHM);
elseif ORG_STRUC.molecule
    ORG_STRUC.toleranceBestHM = 0.10;
end
% maximum of gaussian for antiseed correction
%[nothing, antiSeedsMax] = unix (['./getStuff ' inputFile ' antiSeedsMax 1']);
antiSeedsMax = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsMax -c 1']);
if ~isempty(antiSeedsMax)
    ORG_STRUC.antiSeedsMax = str2num(antiSeedsMax);
end
% sigma of gaussian for antiseed correction
%[nothing, antiSeedsSigma] = unix (['./getStuff ' inputFile ' antiSeedsSigma 1']);
antiSeedsSigma = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsSigma -c 1']);
if ~isempty(antiSeedsSigma)
    ORG_STRUC.antiSeedsSigma = str2num(antiSeedsSigma);
end
% switch the use of order in variation operators (except coormutation) on and off
%[nothing, ordering] = unix (['./getStuff ' inputFile ' ordering_active 1']);
ordering = python_uspex(getPy, ['-f ' inputFile ' -b ordering_active -c 1']);
if ~isempty(ordering)
    ORG_STRUC.ordering = str2num(ordering);
end
% how many generations wait to add extra antiseed if the best structure is not changing
%[nothing, antiSeedsActivation] = unix (['./getStuff ' inputFile ' antiSeedsActivation 1']);
antiSeedsActivation = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsActivation -c 1']);
if ~isempty(antiSeedsActivation)
    ORG_STRUC.antiSeedsActivation = str2num(antiSeedsActivation);
end


if ORG_STRUC.varcomp == 1 | ORG_STRUC.dimension == -4 % only varcomp or proteins
    ORG_STRUC.weight=1;
else
    if ORG_STRUC.molecule == 1
      [type, MtypeLIST, numIons]=GetPOP_MOL(ORG_STRUC.numMols);
    else
      numIons = ORG_STRUC.numIons;
    end
    % weight needed for normalisation of the cosine distance between fingerprints
    L = length(numIons);
    S = 0;
    ORG_STRUC.weight = zeros(L*L,1);
    for i = 1:L
        for j = 1:L
            ind = (i-1)*L+j;
            ORG_STRUC.weight(ind) = (numIons(i)*numIons(j));
            S = S + (numIons(i)*numIons(j));
        end
    end
    
    ORG_STRUC.weight = ORG_STRUC.weight/S;
end
