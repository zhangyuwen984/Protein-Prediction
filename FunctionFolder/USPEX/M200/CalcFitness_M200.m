function fitness = CalcFitness_M200()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
for fit_loop = 1:length(POP_STRUC.POPULATION)
   if ORG_STRUC.optType == 1 % enthalpy
       N_atoms = sum(POP_STRUC.POPULATION(fit_loop).numIons);
       fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end)/N_atoms;
   elseif ORG_STRUC.optType == 7 % bandgap
       fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).gap;
   end
end

fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)
for i = 1 : length(fitness)
  if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
     fitness(i) = 100000;
  end
end
