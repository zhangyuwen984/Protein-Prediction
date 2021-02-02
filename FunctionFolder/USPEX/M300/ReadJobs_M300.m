function ReadJobs_M300()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    % check whether calculations are still performed and read out the output if they are done
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
       Step = POP_STRUC.POPULATION(whichInd).Step;
       disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
       if POP_STRUC.POPULATION(whichInd).JobID
           if ORG_STRUC.platform > 0 | (ORG_STRUC.numParallelCalcs > 1)
              disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
           end
           doneOr = checkStatusC(whichInd);
           if doneOr

              if POP_STRUC.POPULATION(whichInd).JobID == 0.01   % reoptOld = 0, survived structure
                 POP_STRUC.POPULATION(whichInd).Step = length([ORG_STRUC.abinitioCode]) + 1;
              else
                 Error =  Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
                 if Error == 0
                    lattice   =POP_STRUC.POPULATION(whichInd).LATTICE;
                    coordinate=POP_STRUC.POPULATION(whichInd).COORDINATES;
                    chanAList =POP_STRUC.POPULATION(whichInd).chanAList;
                    candidate = reduce_GB(lattice,coordinate, chanAList, ORG_STRUC.thicknessS, ORG_STRUC.numIons);
                    POP_STRUC.POPULATION(whichInd).GB_COORDINATES=candidate;
                 end
              end

              POP_STRUC.POPULATION(whichInd).JobID = 0;

              if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
                 POP_STRUC.POPULATION(whichInd).Done = 1;
                 POP_STRUC.POPULATION(whichInd).ToDo = 0;
                 POP_STRUC.POPULATION(whichInd).Folder=0;

               else
                  if POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode])
                     POP_STRUC.POPULATION(whichInd).Done = 1;
                     POP_STRUC.POPULATION(whichInd).ToDo = 0;
                     POP_STRUC.POPULATION(whichInd).Folder=0;
                     POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                     POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
   
                     numIons = POP_STRUC.POPULATION(whichInd).numIons;
                     [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(lattice, coordinate, numIons, ORG_STRUC.atomType);
                     [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
                     POP_STRUC.POPULATION(whichInd).order =  order;
                     POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
                     POP_STRUC.POPULATION(whichInd).GB_order = reduce_order(order, POP_STRUC.POPULATION(whichInd).chanAList);
                     
                     disp('Relaxation is done.')
                     disp(' ')
                     POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(whichInd, atom_fing);
                     POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);

                     POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                     WriteIndividualOutput_M300(whichInd);
                  else
                     %%%%%%%%%%%%%%%makeGB
                     lattice=POP_STRUC.POPULATION(whichInd).LATTICE;
                     coordinate=POP_STRUC.POPULATION(whichInd).COORDINATES;
                  end
              end
              safesave ('Current_POP.mat', POP_STRUC)
           end
       end
    end
end
