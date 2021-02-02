function EA_311()

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
global POOL_STRUC
global OFF_STRUC
global ANTISEEDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1 % this eternal cycle is needed for non parallel calculations, parallelized one will /break out of it
        
        ReadJobs_311();
        SubmitJobs_311();
        
        if (ORG_STRUC.platform > 0) | (ORG_STRUC.numParallelCalcs > 1)
            if sum([POP_STRUC.POPULATION(:).Done])~= length(POP_STRUC.POPULATION)
                [nothing, nothing] = unix('rm still_running');
                fclose ('all');
                quit
            else
                break;
            end
        else
            if sum([POP_STRUC.POPULATION(:).Done]) == length(POP_STRUC.POPULATION)
                break; % break out of eternal cycle when non parallel calculations finished for given generation
            end
        end
    end
    
    disp('status = Local optimisation finished')
    disp(' ')
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp, [alignLine('Local optimization finished') '\n']);
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', POP_STRUC.generation) ) '\n'] );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    good_and_fit = 0;
    for fit_loop = 1 : length(POP_STRUC.POPULATION)
        if POP_STRUC.POPULATION(fit_loop).Enthalpies(end) < 90000
            good_and_fit = good_and_fit + 1;
        end
    end
    
    if good_and_fit < floor(length(POP_STRUC.POPULATION)/3)
        fprintf(fp,'Too many structures have errors or failed the constraints after optimization.\n');
        fprintf(fp,'Please check the input files. The calculation has to stop.\n');
        fprintf(fp,'Possible reasons: badly tuned optimization parameters or unreasonable contraints.\n');
        quit;
    elseif good_and_fit/length(POP_STRUC.POPULATION) < ORG_STRUC.bestFrac
        fprintf(fp,'Some structures have errors or failed the constraints after optimisation,\n');
        fprintf(fp,'bestFrac parameter will be lowered for this generation to discard such structures.\n');
        fprintf(fp,'Please check the input files. Results may not be reliable. \n');
        fprintf(fp,'Possible reasons: high bestFrac, bad optimization parameters or unreasonable contraints.\n');
    end
    
    if ~((ORG_STRUC.pickedUP == 1) & (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        fitness = CalcFitness_311();
        cor_coefficient = Correlation(fitness);
        fprintf(fp, [alignLine( sprintf('Correlation coefficient = %.4f', cor_coefficient) ) '\n'] );
        
        AntiSeedsCorrection(fitness);
        fitness = FitnessRankingCorrection(fitness);
        
        
        POOL_STRUC.Composition_ranking = zeros(length(POOL_STRUC.Composition_ranking),1);
        POOL_STRUC.POPULATION  = [];
        POOL_STRUC.POPULATION  = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons', {}, 'numBlocks',{}, 'numMols', {}, 'order', {}, 'enthalpy', {}, 'Number', {});
        
        compositionStatistic(ORG_STRUC.firstGeneSplitAll);
        WriteGenerationOutput_311(fitness);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if POP_STRUC.generation >= ORG_STRUC.numGenerations
        Finish();
    else
        
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,'\n');
        if ~isempty(POP_STRUC.convex_hull)
            total_comp = sum(POOL_STRUC.Composition_ranking(:)>0);
            fprintf(fp,'The number of compositions in the current selection pool: %4d\n', total_comp);
            total_pop = length(POOL_STRUC.POPULATION);
            fprintf(fp,'The number of   structures in the current selection pool: %4d\n', total_pop);
            fprintf(fp,'\n');
            fprintf(fp,'Compositions      N_Structures     BestEnthalpies(eV/atom)  fitness    Surviving Gen\n');
            for i=1:length(POOL_STRUC.Composition_ranking)
                comp = num2str(POOL_STRUC.Composition_ratio(i,:), '%6.3f');
                rank = POOL_STRUC.Composition_ranking(i);
                best = POOL_STRUC.Composition_Bestenthalpy(i);
                surv = POOL_STRUC.Composition_surviving(i);
                fit = POOL_STRUC.Composition_fitness(i);
                if rank > 0
                    fprintf(fp,'%-18s    %4d         %10.4f     %16.3f       %4d\n', comp, rank, best, fit, surv);
                end
            end
        end
        
        fprintf(fp,'\n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        POOL2MOL();
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WriteGenerationBackup();
        update_STUFF('INPUT.txt', good_and_fit/length(POP_STRUC.POPULATION), POP_STRUC.ranking);
        poolsize = length(POOL_STRUC.POPULATION);
        ORG_STRUC.tournament = zeros(poolsize,1);
        ORG_STRUC.tournament(poolsize) = 1;
        for loop = 2:poolsize
            ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2) + loop^2;
        end
        
        if POP_STRUC.generation == 1
            if poolsize > ORG_STRUC.initialPopSize
                ORG_STRUC.tournament = zeros(ORG_STRUC.initialPopSize,1);
                ORG_STRUC.tournament(ORG_STRUC.initialPopSize) = 1;
                for loop = 2:ORG_STRUC.initialPopSize
                    ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2)+loop^2;
                end
            end
        end
        
        CreateCalcFolder();
        
        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %just make it has same structure
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);
        
        
        % some generations only softmutated, if needed
        if (ORG_STRUC.softMutOnly(POP_STRUC.generation+1) == 1)
            ORG_STRUC.howManyMutations = 0;
            ORG_STRUC.howManyPermutations = 0;
            ORG_STRUC.howManyRand = 0;
            ORG_STRUC.howManyAtomMutations = length(POP_STRUC.ranking) - POP_STRUC.bad_rank;
            ORG_STRUC.howManyOffsprings = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_311', 'Random_311', 'Permutation_311', ...
            'SoftModeMutation_311', 'LatMutation_311', 'Rotation_311'};
        Num_Opera = [ORG_STRUC.howManyOffsprings, ORG_STRUC.howManyRand, ORG_STRUC.howManyPermutations, ...
            ORG_STRUC.howManyAtomMutations, ORG_STRUC.howManyMutations, ORG_STRUC.howManyRotations];
        count = 0;
        for i = 1 : length(Num_Opera)
            for j = 1:Num_Opera(i)
                count = count + 1;
                eval([Operation{i} '(' num2str(count) ')']);
            end
        end
        
        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        
        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;
        
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Variation operators applied') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,'            %4d structures produced by heredity     \n', ORG_STRUC.howManyOffsprings);
        fprintf(fp,'            %4d structures produced by random       \n', ORG_STRUC.howManyRand);
        fprintf(fp,'            %4d structures produced by softmutation \n', ORG_STRUC.howManyAtomMutations);
        fprintf(fp,'            %4d structures produced by permutation  \n', ORG_STRUC.howManyPermutations);
        fprintf(fp,'            %4d structures produced by latmutation  \n', ORG_STRUC.howManyMutations);
        fprintf(fp,'            %4d structures produced by rotmutation  \n', ORG_STRUC.howManyRotations);
        
        addon_diff = KeepBestStructures();
        fprintf(fp,'            %4d structures kept as best from the previous generation\n', addon_diff);
        
        POP_STRUC = OFF_STRUC;
        num = pick_Seeds();
        fprintf(fp,'            %4d Seeds structures are added from Seeds/POSCARS_%3d \n', num, POP_STRUC.generation);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to the new generation relaxation') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,'            %4d parallel calculations are performed simutaneously\n', ORG_STRUC.numParallelCalcs);
        fprintf(fp, [alignLine('-', 0) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
        fprintf(fp,'  ID   Origin      Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY\n');
        WriteGenerationStart();
        Start_POP_311();
        
        safesave ('Current_POP.mat', POP_STRUC)
        safesave ('Current_ORG.mat', ORG_STRUC)
        disp(' ');
        disp('New generation produced');
        
    end
end

fclose(fp);

