function EA_300()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global ANTISEEDS
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

ToStop = 0;

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    %if ORG_STRUC.fixRndSeed > 0
    %	rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while  ~ORG_STRUC.currentGenDone  % this eternal cycle is needed for non parallel calculations, parallelized one will /break out of it
        
        ReadJobs_300();
        SubmitJobs_300();
        
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
        disp('status = Local optimization finished')
        disp(' ')
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Local optimization finished') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', POP_STRUC.generation) ) '\n'] );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ORG_STRUC.abinitioCode(1) > 0
            good_and_fit = 0;
            for fit_loop = 1 : length(POP_STRUC.POPULATION)
                if POP_STRUC.POPULATION(fit_loop).Enthalpies(end) < 90000
                    good_and_fit = good_and_fit + 1;
                end
            end
        else
            good_and_fit = length(POP_STRUC.POPULATION);
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
            fprintf(fp,'Please check the input files. Results may not be reliable. \n');
            fprintf(fp,'Possible reasons: high bestFrac, bad optimization parameters or unreasonable contraints.\n');
        end
        
        POP_STRUC.good_and_fit = good_and_fit;
        
        if ~((ORG_STRUC.pickedUP == 1) & (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            fitness = CalcFitness_300();
            cor_coefficient = Correlation(fitness);
            fprintf(fp, [alignLine( sprintf('Correlation coefficient = %.4f', cor_coefficient) ) '\n'] );
            
            fitness = AntiSeedsCorrection(fitness);
            fitness = FitnessRankingCorrection(fitness);
            
            compositionStatistic(ORG_STRUC.firstGeneSplitAll);
            
            WriteGenerationOutput_300(fitness);

            if(ORG_STRUC.mlPeerFilter > 0 && ORG_STRUC.mlMinPopulation <= length(USPEX_STRUC.POPULATION))
                mlFilterTrain();
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ORG_STRUC.currentGenDone=1;
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
        ToStop = 0;
        if ~((ORG_STRUC.pickedUP == 1) & (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            [ToStop, Message] = StopRun(fitness, POP_STRUC.generation, ...
                ORG_STRUC.numGenerations, ORG_STRUC.stopCrit, ORG_STRUC.stopFitness);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else
        waitForCOPEX();
        
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,'\n');
        
        WriteGenerationBackup();
        
        update_STUFF('INPUT.txt', POP_STRUC.good_and_fit/length(POP_STRUC.POPULATION), POP_STRUC.ranking);
        
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
        
        % A local flag for performing mlPeerFilter
	do_ml_filtering = false;
        if isfield(ORG_STRUC, 'mlPeerFilter') && (ORG_STRUC.mlPeerFilter > 0)
	  min_generation = ceil(ORG_STRUC.mlMinPopulation/ORG_STRUC.populationSize);
          if(POP_STRUC.generation >= min_generation)
	    do_ml_filtering = true;
	  end
	end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_300', 'Random_300', 'Permutation_300', 'LatticeMutation_300', 'SoftModeMutation_300'};
        Num_Opera = [ORG_STRUC.howManyOffsprings, ORG_STRUC.howManyRand, ORG_STRUC.howManyPermutations, ...
            ORG_STRUC.howManyMutations,  ORG_STRUC.howManyAtomMutations];
        count = 1;
        count_end = ORG_STRUC.populationSize;
        if do_ml_filtering
	  count_end = ORG_STRUC.populationSize*2;
          Num_Opera = Num_Opera*2;
	end

        for i = 1 : length(Num_Opera)
            for j = 1:Num_Opera(i)
                if count <= count_end
                    eval([Operation{i} '(' num2str(count) ')']);
                    count = length(OFF_STRUC.POPULATION) + 1; %SOFTMutation might give two structures
                end
            end
        end
        if exist('spinOperation')
            spinOperation(count-1);   %NEED TO CHECK.........
        end
        
        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        
        if do_ml_filtering
	  OFF_STRUC.POPULATION = mlFilterPredict(OFF_STRUC.POPULATION, ORG_STRUC.populationSize);
	end

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
        fprintf(fp,'            %4d structures produced by spinmutation \n', ORG_STRUC.howManySpinmutations);
        
        addon_diff = KeepBestStructures();
        fprintf(fp,'            %4d structures kept as best from the previous generation\n', addon_diff);
        [addon_copex,OFF_STRUC] = importCOPEXStructures(OFF_STRUC, POP_STRUC.generation);
        fprintf(fp,'            %4d structures imported from the other USPEX Calculations\n', addon_copex);

        POP_STRUC = OFF_STRUC;
        num = pick_Seeds();
        fprintf(fp,'            %4d Seeds structures are added from Seeds/POSCARS\n', num);
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
        Start_POP_300();
        
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
        
    end
end
fclose(fp);
