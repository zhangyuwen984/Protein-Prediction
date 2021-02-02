function [accepted, newInd, Na] = RemoveStruc(comp, fitness, fit_max, ranking, USPEX_STRUC, generation)

% now remove all structures above threshold from consideration
N = length(fitness);
accepted = zeros(1, N);
newInd   = zeros(1, N);
Na = 0; % number of accepted structures
for i = 1 : N
    if fitness(ranking(i)) <= fit_max
        Na = Na + 1;
        if fitness(ranking(i)) >= -0.05 %remove structures with wrong energy
            accepted(ranking(i)) = 1;
            if USPEX_STRUC.POPULATION(ranking(i)).gen == generation
                newInd(ranking(i)) = 1;
            end
        end
    else
        break;
    end
end
% Remove all duplicates:
%--------------------------------------------------------------------------
% First step: fine filtering:
for i = 2 : Na
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j = 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
                    if SameStructure(ranking(i), ranking(j), USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                        break;
                    end
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
% Second step: rough filtering using A_order:
for i = 2 : Na
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j = 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
                    if SameStructure_order(ranking(i), ranking(j), USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                        accepted(ranking(i)) = 0;
                        break;
                    end
                end
            end
        end
    else
        accepted(ranking(i)) = 0;
    end
end
