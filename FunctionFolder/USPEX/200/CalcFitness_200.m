function fitness = CalcFitness_200()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));

for fit_loop = 1:length(POP_STRUC.POPULATION)
     fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end);
     fitness(fit_loop)=fitness(fit_loop)/prod(POP_STRUC.POPULATION(fit_loop).cell);
end


fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)
