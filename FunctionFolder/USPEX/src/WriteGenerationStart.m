function WriteGenerationStart()

global POP_STRUC

POP_STRUC.ranking = 0;
POP_STRUC.bad_rank = 0;

POP_STRUC.current_dir = [POP_STRUC.resFolder '/generation' num2str(POP_STRUC.generation)];
AuxDir = [POP_STRUC.resFolder '/AuxiliaryFiles'];
if ~exist(AuxDir)
    mkdir(AuxDir);
end
if ~exist(POP_STRUC.current_dir)
    mkdir(POP_STRUC.current_dir);
end
