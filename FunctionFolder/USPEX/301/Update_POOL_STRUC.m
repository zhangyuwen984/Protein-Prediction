function count = Update_POOL_STRUC(count, j, enth, fitness, POPULATION)

global POOL_STRUC
global ORG_STRUC

stopCrit = ORG_STRUC.stopCrit;
bestHM   = max(ORG_STRUC.keepBestHM, 4);       %in the case of a very small size

symg      = POPULATION.symg;
num       = POPULATION.numIons;
lat       = POPULATION.LATTICE;
coords    = POPULATION.COORDINATES;
order     = POPULATION.order;
numBlocks = POPULATION.numBlocks;
volume    = POPULATION.Vol/sum(num);
ratio = num/sum(num);

store = 0;
newcomposition=1;
for ind=1:size(POOL_STRUC.Composition_ratio,1)
    if sum(abs(ratio-POOL_STRUC.Composition_ratio(ind,:))) < 0.01  % old composition
        if POOL_STRUC.Composition_surviving(ind)<stopCrit %we don't allow
            E_ranking=POOL_STRUC.Composition_ranking(ind) + 1;  %Composition_ranking, how many structures
            if E_ranking == 1
                enthalpy1 = enth(j);
                POOL_STRUC.Composition_fitness(ind) = fitness(j);
                if enthalpy1 < POOL_STRUC.Composition_Bestenthalpy(ind)
                    POOL_STRUC.Composition_Bestenthalpy(ind) = enthalpy1;   %update bestenthalpy
                    POOL_STRUC.Composition_surviving(ind) = 1;
                else
                    POOL_STRUC.Composition_surviving(ind) = POOL_STRUC.Composition_surviving(ind)+1;
                end
            end
            if E_ranking < bestHM
                store = 1;
                POOL_STRUC.Composition_ranking(ind) = E_ranking;
            end
        else
            POOL_STRUC.Composition_ranking(ind) = -1 ; %means we do not want this composition
        end
        newcomposition = 0;
        break;
    end
end

if newcomposition==1;  %new composition
    POOL_STRUC.Composition_ratio(end+1,:)= ratio;
    POOL_STRUC.Composition_ranking(end+1) = 1;
    POOL_STRUC.Composition_surviving(end+1)= 1;
    POOL_STRUC.Composition_fitness(end+1)= fitness(j);
    POOL_STRUC.Composition_Bestenthalpy(end+1)= enth(j);
    store = 1;
end

if store == 1  % we select only a few structures for each composition
    count     = count + 1;
    POOL_STRUC.POPULATION(count).LATTICE = lat;
    POOL_STRUC.POPULATION(count).COORDINATES = coords;
    POOL_STRUC.POPULATION(count).numIons = num;
    POOL_STRUC.POPULATION(count).numBlocks = numBlocks;
    POOL_STRUC.POPULATION(count).order = order;
    POOL_STRUC.POPULATION(count).enthalpy = enth(j);
    POOL_STRUC.POPULATION(count).Number = j;
end
