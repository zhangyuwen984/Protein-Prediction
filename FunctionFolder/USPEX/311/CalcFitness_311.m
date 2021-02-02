function fitness = CalcFitness_311()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
if ORG_STRUC.optType == 1 % enthalpy
    if size(ORG_STRUC.numIons,1)==1         % two identical blocks
        POP_STRUC.convex_hull = [];
        for fit_loop = 1:length(POP_STRUC.POPULATION)
            fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end)/sum(POP_STRUC.POPULATION(fit_loop).numBlocks);
        end
    else
        fitness = update_convex_hull;
    end
else
    for i = 1:length(POP_STRUC.POPULATION)
        if POP_STRUC.POPULATION(i).Enthalpies(end) < 9999
            if ORG_STRUC.optType == 2 % volume
                fitness(i) = det(POP_STRUC.POPULATION(i).LATTICE)/sum(POP_STRUC.POPULATION(i).numIons);
            elseif ORG_STRUC.optType == 3 % hardness
                fitness(i) = -1*POP_STRUC.POPULATION(i).hardness;
            elseif ORG_STRUC.optType == 4 % Structural Order
                fitness(i) = -1*POP_STRUC.POPULATION(i).S_order;
            elseif ORG_STRUC.optType == 6 % dielectric tensor (or rather Tr()/3)
                fitness(i) = -1*sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
            elseif ORG_STRUC.optType == 7 % bandgap
                fitness(i) = -1*POP_STRUC.POPULATION(i).gap;
            elseif ORG_STRUC.optType == 8 % Tr(dielectric tensor)/3 multiplied by gap^2 value
                Egc = 4; %critical bandgap for semiconductor and insulator, eV
                if POP_STRUC.POPULATION(i).gap >= Egc
                    fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^2; %Insulator
                else
                    fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^6; %Semiconductor
                end
            elseif ORG_STRUC.optType == 9
                fitness(i) = -1*POP_STRUC.POPULATION(i).mag_moment;
            elseif ORG_STRUC.optType == 10 % structure entropy
                fitness(i) = -1*POP_STRUC.struc_entr(i);
            end
        end
    end
end


fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)
for i = 1 : length(fitness)
    if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
        fitness(i) = 100000;
    end
end
