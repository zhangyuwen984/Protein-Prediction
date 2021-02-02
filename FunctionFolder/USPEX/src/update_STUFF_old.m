function update_STUFF_old(inputFile, goodFrac, ranking)
% $Rev: 1006 $
% $Author: qian $
% $Date: 2015-04-25 15:36:13 +0400 (Sat, 25 Apr 2015) $

global ORG_STRUC

% sense. So in that case we stick to pure heredity.
if size(ORG_STRUC.numIons,2)==1
    ORG_STRUC.fracPerm = 0;
    ORG_STRUC.fracGene = ORG_STRUC.fracGene+ORG_STRUC.fracPerm;
end

howManyOffsprings     = round(ORG_STRUC.populationSize*ORG_STRUC.fracGene);
howManyPermutations   = round(ORG_STRUC.populationSize*ORG_STRUC.fracPerm);
howManyAtomMutations  = round(ORG_STRUC.populationSize*ORG_STRUC.fracAtomsMut);
howManyRand           = round(ORG_STRUC.populationSize*ORG_STRUC.fracRand);
howManyRotations      = round(ORG_STRUC.populationSize*ORG_STRUC.fracRotMut);      % molecules & proteins
howManyTransmutations = round(ORG_STRUC.populationSize*ORG_STRUC.fracTrans);       % varcom only
howManySpinmutations  = round(ORG_STRUC.populationSize*ORG_STRUC.fracSpin);       % spin operation
howManySecSwitch      = round(ORG_STRUC.populationSize*ORG_STRUC.fracSecSwitch);   % proteins only
howManyShiftBorder    = round(ORG_STRUC.populationSize*ORG_STRUC.fracShiftBorder); % proteins only

sum_offsprings = howManyOffsprings     + howManyPermutations  + ...
                 howManyTransmutations + howManyAtomMutations + ...
                 howManyRand           + howManyRotations     + ...
                 howManySpinmutations  + ...
                 howManySecSwitch      + howManyShiftBorder;

if sum_offsprings > ORG_STRUC.populationSize
    howManyOffsprings     = round((howManyOffsprings    *ORG_STRUC.populationSize)/sum_offsprings);
    howManyPermutations   = round((howManyPermutations  *ORG_STRUC.populationSize)/sum_offsprings);
    howManyAtomMutations  = round((howManyAtomMutations *ORG_STRUC.populationSize)/sum_offsprings);
    howManyRand           = round((howManyRand          *ORG_STRUC.populationSize)/sum_offsprings);
    howManyTransmutations = round((howManyTransmutations*ORG_STRUC.populationSize)/sum_offsprings);
    howManyRotations      = round((howManyRotations     *ORG_STRUC.populationSize)/sum_offsprings);
    howManySpinmutations  = round((howManySpinmutations *ORG_STRUC.populationSize)/sum_offsprings);
    howManySecSwitch      = round((howManySecSwitch     *ORG_STRUC.populationSize)/sum_offsprings);
    howManyShiftBorder    = round((howManyShiftBorder   *ORG_STRUC.populationSize)/sum_offsprings);
    sum_offsprings        = howManyOffsprings     + howManyPermutations  + ...
                            howManyTransmutations + howManyAtomMutations + ...
                            howManyRand           + howManyRotations     + ...
                            howManySpinmutations  + ...
                            howManySecSwitch      + howManyShiftBorder;
end

howManyleft = max([ORG_STRUC.populationSize-sum_offsprings, 0]);
if ~ORG_STRUC.constLattice
    howManyMutations = howManyleft;  %lat_Mutation
else
    howManyOffsprings = howManyOffsprings + howManyleft;
    howManyMutations = 0;  %lat_Mutation
end

ORG_STRUC.howManyMutations     = howManyMutations;
ORG_STRUC.howManyPermutations  = howManyPermutations;
ORG_STRUC.howManyAtomMutations = howManyAtomMutations;
ORG_STRUC.howManyRand          = howManyRand;
ORG_STRUC.howManyRotations     = howManyRotations;      % molecules & proteins
ORG_STRUC.howManyOffsprings    = howManyOffsprings;
ORG_STRUC.howManyTrans         = howManyTransmutations;
ORG_STRUC.howManySpinmutations = howManySpinmutations;% spin
ORG_STRUC.howManySecSwitch     = howManySecSwitch;      % proteins
ORG_STRUC.howManyShiftBorder   = howManyShiftBorder;    % proteins

end
