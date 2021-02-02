function fitness = CalcFitness_M400()
% $Rev: 696 $
% $Author: maxim $
% $Date: 2014-10-30 21:07:20 +0400 (Thu, 30 Oct 2014) $

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
%Enthalpy
for fit_loop = 1:length(POP_STRUC.POPULATION)
    fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end);
end

fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)
for i = 1 : length(fitness)
    if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
        fitness(i) = 100000;
    end
end
