function correlation_coefficient = Correlation(fitness)

global ORG_STRUC
global POP_STRUC

% correlation between fitness and order:
ORG_STRUC.cor_dir = 1;
ORG_STRUC.correlation_coefficient = 0;
j = 1;
FvsO = zeros(1,2);
for i = 1 : length(POP_STRUC.POPULATION)
    if (fitness(i) < 99999) && (POP_STRUC.DoneOrder(i) > 0)
        a_o = sum(POP_STRUC.POPULATION(i).order)/sum(POP_STRUC.POPULATION(i).numIons);
        if j == 1
            FvsO(1,1) = a_o;
            FvsO(1,2) = fitness(i);
        else
            FvsO = vertcat(FvsO, [a_o fitness(i)]);
        end
        j = j + 1;
    end
end
if exist('corrcoef')
    corr_tmp = corrcoef(FvsO);
else
    corr_tmp = 0;
end
if length(corr_tmp)>1
    correlation_coefficient = corr_tmp(1,2);
else
    correlation_coefficient = 0;
end

ORG_STRUC.correlation_coefficient = correlation_coefficient;
ORG_STRUC.cor_dir = sign(correlation_coefficient); % > 0 - correlation, < 0 - anticorrelation
