function update_STUFF_M400(inputFile, goodFrac, ranking)
% $Rev: 1115 $
% $Author: mrakitin $
% $Date: 2015-08-18 01:44:25 +0400 (Tue, 18 Aug 2015) $

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Some parameters can be updated from INPUT.txt %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
createORG_Fingerprint(inputFile);
createORG_EA(inputFile);
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%[nothing, numParallelCalcs] = unix (['./getStuff ' inputFile ' numParallelCalcs 1']);
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
    ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fraction of each operator can be updated as below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

howManyProliferate = min( [round((length(POP_STRUC.ranking))*ORG_STRUC.bestFrac), length(POP_STRUC.ranking)-POP_STRUC.bad_rank ] );
if ~ORG_STRUC.AutoFrac
    update_STUFF_old(inputFile, goodFrac, ranking);
else
    N = zeros(3,5); % rand, here, rot, ssw, fbord      % it was (3,7)
    % parents, current pop, new pop
    frand = 0; fhere = 0; frot = 0; fssw = 0; fshb = 0;
    gen = POP_STRUC.generation;
    if gen > 1
        prev1_gen  = USPEX_STRUC.GENERATION(gen-1).ID;  %
        
        if gen == 2
            prev2_gen  = 0;
        else
            prev2_gen  = USPEX_STRUC.GENERATION(gen-2).ID;  %
        end
        for i = 1:howManyProliferate %for all the good Structures
            C_ID = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).Number;
            accept = 1;
            for j = prev2_gen + 1 : prev1_gen
                if SameStructure_order(j, C_ID, USPEX_STRUC)  %remove the structures which are same from previous gen
                    accept = 0;
                    break;
                end
            end
            
            if accept
                tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
                if ~isempty(findstr(tmp, 'Random'))
                    N(1,1) = N(1,1) + 1;
                elseif ~isempty(findstr(tmp, 'Heredity'))
                    N(1,2) = N(1,2) + 1;
                elseif ~isempty(findstr(tmp, 'Rotation'))
                    N(1,3) = N(1,3) + 1;
                elseif ~isempty(findstr(tmp, 'SecSwitch'))
                    N(1,4) = N(1,4) + 1;
                elseif ~isempty(findstr(tmp, 'ShiftBorder'))
                    N(1,5) = N(1,5) + 1;
                end
            end
        end
        
        for i = 1:length(POP_STRUC.POPULATION) %for all the good Structures
            tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
            if ~isempty(findstr(tmp, 'Random'))
                N(2,1) = N(2,1) + 1;
            elseif ~isempty(findstr(tmp, 'Heredity'))
                N(2,2) = N(2,2) + 1;
            elseif ~isempty(findstr(tmp, 'Rotation'))
                N(2,3) = N(2,3) + 1;
            elseif ~isempty(findstr(tmp, 'SecSwitch'))
                N(2,4) = N(2,4) + 1;
            elseif ~isempty(findstr(tmp, 'ShiftBorder'))
                N(2,5) = N(2,5) + 1;
            end
        end
        
        count = 0;
        sum1 = 0;
        for i = 1:5     % it was 7
            if N(2,i) > 0
                X(i) = N(1,i)/N(2,i);
                count = count + 1;
                sum1 = sum1 + X(i);
            end
        end
        
        X_mean = sum1/count;
        for i = 1:5     % it was 7
            if N(2,i) >0
                N(3,i) = round(N(1,i)*(X(i)/X_mean+1)/2);
            end
        end
        
        T_N = sum(N(2,:));
        Ar  = N(2,1)/(T_N*ORG_STRUC.fracRand);
        N(3,1) = N(3,1)/Ar;
        T = sum(N(3,:));
        frand = N(3,1)/T;
        fhere = N(3,2)/T;
        frot  = N(3,3)/T;
        fssw  = N(3,4)/T;
        fshb  = N(3,5)/T;
    end    % if gen > 0
    
    frand = max((ORG_STRUC.fracRand        + frand)/2, 0.10);
    fhere = max((ORG_STRUC.fracGene        + fhere)/2, 0.10);
    frot  = max((ORG_STRUC.fracRotMut      + frot )/2, 0.10);
    fssw  = max((ORG_STRUC.fracSecSwitch   + fssw )/2, 0.10);
    fshb  = max((ORG_STRUC.fracShiftBorder + fshb )/2, 0.10);
    
    T_frac = frand + fhere + frot + fssw + fshb;
    ORG_STRUC.fracRand        = frand/T_frac;
    ORG_STRUC.fracGene        = fhere/T_frac;
    ORG_STRUC.fracRotMut      = frot /T_frac;
    ORG_STRUC.fracSecSwitch   = fssw /T_frac;
    ORG_STRUC.fracShiftBorder = fshb /T_frac;
    
    ORG_STRUC.howManyRand         = round(ORG_STRUC.populationSize*ORG_STRUC.fracRand);
    ORG_STRUC.howManyOffsprings   = round(ORG_STRUC.populationSize*ORG_STRUC.fracGene);
    ORG_STRUC.howManyRotations    = round(ORG_STRUC.populationSize*ORG_STRUC.fracRotMut);
    ORG_STRUC.howManySecSwitch    = round(ORG_STRUC.populationSize*ORG_STRUC.fracSecSwitch);
    ORG_STRUC.howManyShiftBorder  = round(ORG_STRUC.populationSize*ORG_STRUC.fracShiftBorder);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% The tournement routine is defined here%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THIS IS IMPORTANT and subject to change. tournament specifies the
% probability of each proliferating individual to be selected as parent or
% mutation template
% the following ORG_STRUC.tournament is quadratic
ORG_STRUC.tournament = zeros(howManyProliferate,1);
ORG_STRUC.tournament(howManyProliferate) = 1;
for loop = 2:howManyProliferate
    ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2) + loop^2;
end

if POP_STRUC.generation == 1
    if howManyProliferate > ORG_STRUC.initialPopSize
        if ORG_STRUC.initialPopSize==0
            ORG_STRUC.initialPopSize=1;
        end
        ORG_STRUC.tournament = zeros(ORG_STRUC.initialPopSize,1);
        ORG_STRUC.tournament(ORG_STRUC.initialPopSize) = 1;
        for loop = 2:ORG_STRUC.initialPopSize
            ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2)+loop^2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%% END defining tournement routine %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(ORG_STRUC.softMutOnly,1)<ORG_STRUC.numGenerations+2
    ORG_STRUC.softMutOnly = zeros(ORG_STRUC.numGenerations + 1,1);
end

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);

fprintf(fp,'The fittest %2d percent of the population used to produce next generation\n', ORG_STRUC.bestFrac*100);
fprintf(fp,'    fraction of generation produced by heredity        :     %4.2f\n', ORG_STRUC.fracGene);
fprintf(fp,'    fraction of generation produced by random          :     %4.2f\n', ORG_STRUC.fracRand);
fprintf(fp,'    fraction of generation produced by rotmutation     :     %4.2f\n', ORG_STRUC.fracRotMut);
fprintf(fp,'    fraction of generation produced by secondary switch:     %4.2f\n', ORG_STRUC.fracSecSwitch);
fprintf(fp,'    fraction of generation produced by shift border    :     %4.2f\n', ORG_STRUC.fracShiftBorder);

fclose(fp);

% this is average enthalpy of all selected structures
ORG_STRUC.averageEnergy = 0;
for i = 1 : howManyProliferate
    ORG_STRUC.averageEnergy = ORG_STRUC.averageEnergy + POP_STRUC.POPULATION(ranking(i)).Enthalpies(end);
end
ORG_STRUC.averageEnergy = ORG_STRUC.averageEnergy/howManyProliferate;

safesave ('Current_ORG.mat', ORG_STRUC)
