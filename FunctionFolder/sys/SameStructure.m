function same = SameStructure(ID1, ID2, POP)
% To check if two structures are identical
% weight:  for fingerprint analysis
% all the structure are loaded from USPEX_STRUC
% Last updated by Qiang Zhu (2014/02/18)

global ORG_STRUC

numIons1 = POP.POPULATION(ID1).numIons;
numIons2 = POP.POPULATION(ID2).numIons;

if sameComposition(numIons1, numIons2)
    weight     = ORG_STRUC.weight;
    toleranceF = ORG_STRUC.toleranceFing;
    if ORG_STRUC.molecule
       toleranceE = 0.005;
       toleranceO = 0.002;
    else
       toleranceE = 0.015;
       toleranceO = 0.075;
    end
    
    E1 = POP.POPULATION(ID1).Enthalpies(end)/sum(numIons1);
    E2 = POP.POPULATION(ID2).Enthalpies(end)/sum(numIons2);
    f1 = POP.POPULATION(ID1).FINGERPRINT;
    f2 = POP.POPULATION(ID2).FINGERPRINT;
    ao1 = mean(POP.POPULATION(ID1).order);
    ao2 = mean(POP.POPULATION(ID2).order);

    dist_F = cosineDistance(f1, f2, weight);
    dist_O = abs(ao1-ao2);
    diff   = abs(E1-E2);
    if (ORG_STRUC.dimension == 0) || (ORG_STRUC.dimension == 2)
        if (dist_F < toleranceF) && (diff < toleranceE)
            same = 1;
        else
            same = 0;
        end
    else
        if ((dist_F < toleranceF) || (dist_O < toleranceO))...
           && (diff < toleranceE)
            same = 1;
        else
            same = 0;
        end
    end
else
    same = 0;
end
