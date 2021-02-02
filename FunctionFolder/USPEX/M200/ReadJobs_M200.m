function ReadJobs_M200()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    % check whether calculations are still performed and read out the output if they are done
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
       Step = POP_STRUC.POPULATION(whichInd).Step;
       disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
       if POP_STRUC.POPULATION(whichInd).JobID
           if (ORG_STRUC.platform > 0) | (ORG_STRUC.numParallelCalcs > 1)
              disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
           end
           doneOr = checkStatusC(whichInd);
           if doneOr
              if POP_STRUC.POPULATION(whichInd).JobID == 0.01   % reoptOld = 0, survived structure
                 POP_STRUC.POPULATION(whichInd).Step = length([ORG_STRUC.abinitioCode]) + 1;
              elseif POP_STRUC.POPULATION(whichInd).JobID ~= 0.02   % no optimization
                 Error = Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
                 if Error == 0
                    lattice=POP_STRUC.POPULATION(whichInd).LATTICE;
                    coordinate=POP_STRUC.POPULATION(whichInd).COORDINATES;
                    [lat,candidate]=reduce2D(lattice,coordinate,ORG_STRUC.thicknessS+3);
                    POP_STRUC.POPULATION(whichInd).LATTICE_2D=lat;
                    POP_STRUC.POPULATION(whichInd).COORDINATES_2D=candidate;
                    POP_STRUC.POPULATION(whichInd).Volume_2D=det(lat);
                    Step = POP_STRUC.POPULATION(whichInd).Step;
                    if Step < length ([ORG_STRUC.abinitioCode])  %for the final step, we don't make new!!!
                       [lattice, coordinate] = make2D(lat, candidate, ORG_STRUC.vacuumSize(Step)-3);
                       POP_STRUC.POPULATION(whichInd).LATTICE = lattice;
                       POP_STRUC.POPULATION(whichInd).COORDINATES = coordinate;
                    else
                       if max(coordinate(:,3))-min(coordinate(:,3)) > (1-2/lattice(3,3)) %too thick
                          POP_STRUC.POPULATION(whichInd).Enthalpies(Step-1)=100000;
                       end
                    end
                 end
              end

              POP_STRUC.POPULATION(whichInd).JobID = 0;

              if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
                 POP_STRUC.POPULATION(whichInd).Done = 1;
                 POP_STRUC.POPULATION(whichInd).ToDo = 0;
                 POP_STRUC.POPULATION(whichInd).Folder=0;

               elseif POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode]) %Step is different
                  POP_STRUC.POPULATION(whichInd).Done = 1;
                  POP_STRUC.POPULATION(whichInd).ToDo = 0;
                  POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                  POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
                  LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE_2D;
                  numIons = POP_STRUC.POPULATION(whichInd).numIons;
                  COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES_2D;
                  atomType = ORG_STRUC.atomType;
                  [Ni, V, dist_matrix, typ_i, typ_j, ho, ht] = makeMatrices_2D(LATTICE, COORDINATES, numIons, atomType);
                  [order, FINGERPRINT, atom_fing] = fingerprint_calc_2D(Ni, V, dist_matrix, typ_i, typ_j, numIons, ho, ht);
                  POP_STRUC.POPULATION(whichInd).order =  order;
                  POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
                  POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(whichInd, atom_fing);
                  POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);

                  disp('Relaxation is done.')
                  disp(' ')
                  POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                  WriteIndividualOutput_M200(whichInd);
                  POP_STRUC.POPULATION(whichInd).Folder=0;
              end
              safesave ('Current_POP.mat', POP_STRUC)
           end
       end
    end
end
