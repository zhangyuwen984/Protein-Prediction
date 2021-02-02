function EA_301()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global POOL_STRUC
global ANTISEEDS
global USPEX_STRUC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ~ORG_STRUC.currentGenDone  % this eternal cycle is needed for non parallel calculations, parallelized one will /break out of it
        
        ReadJobs_301();
        SubmitJobs_301();
        
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
    
    if ORG_STRUC.currentGenDone == 0
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
            fprintf(fp,'selection of these structures as parents is possible.\n');
            fprintf(fp,'bestFrac parameter will be lowered for this generation to discard such structures.\n');
            fprintf(fp,'Possible reasons: high bestFrac, bad optimization parameters or unreasonable contraints.\n');
        end
        
        POP_STRUC.good_and_fit = good_and_fit;
        if ~((ORG_STRUC.pickedUP == 1) & (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            fitness = CalcFitness_301();
            cor_coefficient = Correlation(fitness);
            fprintf(fp, [alignLine( sprintf('Correlation coefficient = %.4f', cor_coefficient) ) '\n'] );
            
            fitness = AntiSeedsCorrection(fitness);
            fitness = FitnessRankingCorrection(fitness);
            POOL_STRUC.Composition_ranking = zeros(length(POOL_STRUC.Composition_ranking),1);
            POOL_STRUC.POPULATION  = [];
            POOL_STRUC.POPULATION  = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons', {}, 'numBlocks',{}, 'order', {}, 'enthalpy', {}, 'Number', {});
            
            compositionStatistic(ORG_STRUC.firstGeneSplitAll);
            WriteGenerationOutput_301(fitness);
        end
        
        WriteGenerationBackup();
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 4:  Update the USPEX.mat and POOL.mat with AntiSeeds     %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reupdateStructures_301();
        
        %Step 5: save and quit
        
        ORG_STRUC.currentGenDone=1;
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
        
        if POP_STRUC.generation >= ORG_STRUC.numGenerations
            Finish();
            fclose('all');
            return;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 6: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Start the new generation %%
    waitForCOPEX();
    
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp,'\n');
    if ~isempty(POP_STRUC.convex_hull)
        total_comp = sum(POOL_STRUC.Composition_ranking(:)>0);
        total_pop = length(POOL_STRUC.POPULATION);
        
        fprintf(fp,'The number of compositions in the current selection pool: %4d\n', total_comp);
        fprintf(fp,'The number of   structures in the current selection pool: %4d\n', total_pop);
        fprintf(fp,'\n');
        fprintf(fp,'Compositions      N_Structures     BestEnthalpies(eV/atom)  fitness    Surviving Gen\n');
        for i=1:length(POOL_STRUC.Composition_ranking)
            comp = num2str(POOL_STRUC.Composition_ratio(i,:), '%6.3f');
            rank = POOL_STRUC.Composition_ranking(i);
            best = POOL_STRUC.Composition_Bestenthalpy(i);
            surv = POOL_STRUC.Composition_surviving(i);
            fit  = POOL_STRUC.Composition_fitness(i);
            if rank > 0
                fprintf(fp,'%-18s    %4d         %10.4f     %16.3f       %4d\n', comp, rank, best, fit, surv);
            end
        end
    end
    
    fprintf(fp,'\n');
    
    update_STUFF('INPUT.txt', POP_STRUC.good_and_fit/length(POP_STRUC.POPULATION), POP_STRUC.ranking);
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
    if POP_STRUC.generation < ORG_STRUC.numGenerations
        if (ORG_STRUC.softMutOnly(POP_STRUC.generation+1) == 1)
            ORG_STRUC.howManyMutations = 0;
            ORG_STRUC.howManyPermutations = 0;
            ORG_STRUC.howManyRand = 0;
            ORG_STRUC.howManyAtomMutations = length(POP_STRUC.ranking) - POP_STRUC.bad_rank;
            ORG_STRUC.howManyOffsprings = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Operation = {'Heredity_301', 'Random_301', 'Permutation_301', 'Transmutation_301', 'SoftModeMutation_301', 'LatticeMutation_301'};
    Num_Opera = [ORG_STRUC.howManyOffsprings, ORG_STRUC.howManyRand, ORG_STRUC.howManyPermutations, ...
        ORG_STRUC.howManyTrans, ORG_STRUC.howManyAtomMutations, ORG_STRUC.howManyMutations];
    count = 0;
    for i = 1 : length(Num_Opera)
        for j = 1:Num_Opera(i)
            count = count + 1;
            eval([Operation{i} '(' num2str(count) ')']);
        end
    end

    if exist('spinOperation')
        spinOperation(count);
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
    fprintf(fp,'            %4d structures produced by transmutation\n', ORG_STRUC.howManyTrans);
    fprintf(fp,'            %4d structures produced by spinmutation \n', ORG_STRUC.howManySpinmutations);
    
    addon_diff = KeepBestStructures();
    fprintf(fp,'            %4d structures kept as best from the previous generation\n',  addon_diff);
    [addon_copex,OFF_STRUC] = importCOPEXStructures(OFF_STRUC, POP_STRUC.generation);
    fprintf(fp,'            %4d structures imported from the other USPEX Calculations\n', addon_copex);
    
    
    POP_STRUC = OFF_STRUC;
    num = pick_Seeds();
    fprintf(fp,'            %4d Seeds structures are added from Seeds/POSCARS\n', num);
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp, [alignLine('Proceeding to the new generation relaxation') '\n']);
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp,'            %4d parallel calculations are performed simultaneously\n', ORG_STRUC.numParallelCalcs);
    fprintf(fp, [alignLine('-', 0) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% THE NEW POPULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Start_POP_301();
    WriteGenerationStart();
    fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
    fprintf(fp,'  ID   Origin     Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY\n');
    
    
    ORG_STRUC.currentGenDone = 0;
    if ( ORG_STRUC.pluginType > 0 ) & ( mod( POP_STRUC.generation, ORG_STRUC.pluginType)==0 )
        ORG_STRUC.startNextGen = 0;
    else
        ORG_STRUC.startNextGen = 1;
    end
    
    safesave ('Current_POP.mat', POP_STRUC)
    safesave ('Current_ORG.mat', ORG_STRUC)
    
    disp(' ');
    disp('New generation produced');
    %Start the new generation
end     %while generation cycle

fclose(fp);

