function SubmitJobs_301()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    DO_NOW = 0;
    if isempty (whichInd)
        stillLeft = find([POP_STRUC.POPULATION(:).ToDo]);
        if ~isempty(stillLeft)
           DO_NOW = stillLeft(1);
           POP_STRUC.POPULATION(DO_NOW).Folder = indic;
           POP_STRUC.POPULATION(DO_NOW).ToDo = 0;
        end
    elseif POP_STRUC.POPULATION(whichInd).JobID == 0
        DO_NOW = whichInd;
    end
    if DO_NOW
        Step = POP_STRUC.POPULATION(DO_NOW).Step;
        if (Step > 1) && (Step < length(ORG_STRUC.abinitioCode)-1)
	    POP_STRUC.POPULATION(DO_NOW).LATTICE = perturbCell(POP_STRUC.POPULATION(DO_NOW).LATTICE,1);
            POP_STRUC.POPULATION(DO_NOW).COORDINATES = perturbCoords(POP_STRUC.POPULATION(DO_NOW).COORDINATES, POP_STRUC.POPULATION(DO_NOW).LATTICE,1); 
        end

        if Step > length ([ORG_STRUC.abinitioCode]) % structures from Best
          POP_STRUC.POPULATION(DO_NOW).JobID = 0.01;
        elseif ORG_STRUC.abinitioCode(Step) == 0   % no optimization at all! (used in order optimization)
          POP_STRUC.POPULATION(DO_NOW).JobID = 0.02;
        else
          cd (['CalcFold' num2str(indic)])
          Write_AbinitCode(ORG_STRUC.abinitioCode(Step), DO_NOW);
          POP_STRUC.POPULATION(DO_NOW).JobID = submitJob(DO_NOW);
          cd ..
        end
          safesave ('Current_POP.mat', POP_STRUC)
    end
end
