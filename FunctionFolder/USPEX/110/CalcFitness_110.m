function fitness = CalcFitness_110()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
if ORG_STRUC.optType == 1 % enthalpy
   for fit_loop = 1:length(POP_STRUC.POPULATION)
       fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end);
   end
end

