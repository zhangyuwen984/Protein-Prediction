function fitness = CalcFitness_201()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
if length(ORG_STRUC.numIons)>1
   fitness = update_convex_hull_surface(); %only 2D
else
   for fit_loop = 1:length(POP_STRUC.POPULATION)
       fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end);
       fitness(fit_loop)=fitness(fit_loop)/sum(POP_STRUC.POPULATION(fit_loop).numIons);
   end
end
