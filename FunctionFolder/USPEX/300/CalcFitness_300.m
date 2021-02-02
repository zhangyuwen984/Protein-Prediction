function fitness = CalcFitness_300()

global POP_STRUC
global ORG_STRUC
fitness = zeros(1,length(POP_STRUC.POPULATION));
if ORG_STRUC.optType == 1 % enthalpy
    for fit_loop = 1:length(POP_STRUC.POPULATION)
        factor = POP_STRUC.POPULATION(fit_loop).numIons/ORG_STRUC.numIons;
        fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end)/factor ;
    end
elseif ORG_STRUC.optType == 2 % volume
    for fit_loop = 1:length(POP_STRUC.POPULATION)
        fitness(fit_loop) = det(POP_STRUC.POPULATION(fit_loop).LATTICE)/sum(POP_STRUC.POPULATION(fit_loop).numIons);
    end
elseif ORG_STRUC.optType == 3 % hardness
    for i = 1 : length(POP_STRUC.POPULATION)
        fitness(i) = -1*POP_STRUC.POPULATION(i).hardness;
    end
elseif ORG_STRUC.optType == 4 % Structural Order
    for fit_loop = 1:length(POP_STRUC.POPULATION)
        fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).S_order;
    end
elseif ORG_STRUC.optType == 5 % average distance to other structures
    for i = 1 : length(POP_STRUC.POPULATION)-1
        for j = i+1 : length(POP_STRUC.POPULATION)
            dist_ij = cosineDistance(POP_STRUC.POPULATION(i).FINGERPRINT, POP_STRUC.POPULATION(j).FINGERPRINT, ORG_STRUC.weight);
            fitness(i) = fitness(i) + dist_ij^2;
            fitness(j) = fitness(j) + dist_ij^2;
        end
    end
    fitness = -sqrt(fitness);
elseif ORG_STRUC.optType == 6 % dielectric tensor (or rather Tr()/3)
    for i = 1 : length(POP_STRUC.POPULATION)
        fitness(i) = -1*sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
    end
elseif ORG_STRUC.optType == 7 % bandgap
    for i = 1 : length(POP_STRUC.POPULATION)
        fitness(i) = -1*POP_STRUC.POPULATION(i).gap;
    end
elseif ORG_STRUC.optType == 8 % Tr(dielectric tensor)/3 multiplied by gap^2 value
    Egc = 4; %critical bandgap for semiconductor and insulator, eV
    for i = 1 : length(POP_STRUC.POPULATION)
        if POP_STRUC.POPULATION(i).gap >= Egc
            fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^2; %Insulator
        else
            fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^6; %Semiconductor
        end
    end
elseif ORG_STRUC.optType == 9
    for fit_loop = 1 : length(POP_STRUC.POPULATION)
        fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).mag_moment;
    end
elseif ORG_STRUC.optType == 10 % structure entropy
    for fit_loop = 1 : length(POP_STRUC.POPULATION)
        factor = POP_STRUC.POPULATION(fit_loop).numIons/ORG_STRUC.numIons;
        fitness(fit_loop) = -1*POP_STRUC.struc_entr(fit_loop)/factor;
    end
elseif ORG_STRUC.optType == 14 % power factor
    for fit_loop = 1 : length(POP_STRUC.POPULATION)
        fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).powerfactor;
    end
elseif (ORG_STRUC.optType > 1100) & (ORG_STRUC.optType < 1112)  % structure elasticProperties
    whichPara= mod(ORG_STRUC.optType,110);
    for i = 1 : length(POP_STRUC.POPULATION)
        if isempty(POP_STRUC.POPULATION(i).elasticProperties) | (POP_STRUC.POPULATION(i).elasticProperties(end)==0)
            fitness(i)=NaN;
        else
            %disp(['Individual ', num2str(i) ]);
            %POP_STRUC.POPULATION(i).elasticProperties
            fitness(i) = -1*POP_STRUC.POPULATION(i).elasticProperties(whichPara);
        end
    end
end

fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)

if ORG_STRUC.abinitioCode(1) > 0
    for i = 1 : length(fitness)
        if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
            fitness(i) = 100000;
        end
    end
end

