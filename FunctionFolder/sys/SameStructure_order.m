function same = SameStructure_order(ID1, ID2, POP)
% To check if two structures are identical
% all the structure are loaded from POP

global ORG_STRUC

numIons1 = POP.POPULATION(ID1).numIons;
numIons2 = POP.POPULATION(ID2).numIons;

if sameComposition(numIons1, numIons2)
    toleranceO = 0.075;
    toleranceE = 0.015;
    
    E1 = POP.POPULATION(ID1).Enthalpies(end)/sum(numIons1);
    E2 = POP.POPULATION(ID2).Enthalpies(end)/sum(numIons2);
    % A_order:
    if ORG_STRUC.molecule == 0
        ao1 = mean(POP.POPULATION(ID1).order);
        ao2 = mean(POP.POPULATION(ID2).order);
    else
        ao1_list = [];
        for i=1:size(POP.POPULATION(ID1).MOLECULES, 2)
            ao1_list = [ao1_list POP.POPULATION(ID1).MOLECULES(i).order];
        end
        ao1 = mean(ao1_list);

        ao2_list = [];
        for i=1:size(POP.POPULATION(ID2).MOLECULES, 2)
            ao2_list = [ao2_list POP.POPULATION(ID2).MOLECULES(i).order];
        end
        ao2 = mean(ao2_list);
    end
    
    dist = abs(ao1 - ao2);
    diff = abs(E1 - E2);
    if (dist < toleranceO) && (diff < toleranceE)
        same = 1;
    else
        same = 0;
    end
else
    same = 0;
end
