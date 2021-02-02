function Random_M400(Ind_No)
% $Rev: 690 $
% $Author: maxim $
% $Date: 2014-10-28 07:25:19 +0400 (Tue, 28 Oct 2014) $

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

angles_num         = size(POP_STRUC.backbone_atoms, 2);

% Fill all necessary variables:
[OFF_STRUC.POPULATION(Ind_No).ANGLES, OFF_STRUC.POPULATION(Ind_No).template_random]   = Random_Init_M400(angles_num);
OFF_STRUC.POPULATION(Ind_No).RESIDUES = POP_STRUC.POPULATION(1).RESIDUES;
OFF_STRUC.POPULATION(Ind_No).Parents = [];
OFF_STRUC.POPULATION(Ind_No).numIons = ORG_STRUC.numIons;
OFF_STRUC.POPULATION(Ind_No).howCome = 'Random';

disp(['Structure ' num2str(Ind_No) ' generated randomly']);

end
