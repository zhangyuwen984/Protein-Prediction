function Write_Tinker(Ind_No)

% Revision: $Rev$
% Author  : $Author$
% Date    : $Date$

global POP_STRUC
global ORG_STRUC
step = POP_STRUC.POPULATION(Ind_No).Step;

if (ORG_STRUC.molecule==1) & (ORG_STRUC.dimension==3)
    try
        [nothing, nothing] = unix(['cat toptions.key_' num2str(step) ' > input.key']);
    catch
        disp('OOUPS, no toptions.key_*');
    end
    writeCalcFilesTINKER_MOL(Ind_No);
else %protein
    if POP_STRUC.generation >= 1
        make_template = strcat('input.make0_', num2str(step));
        % Inserting new coordinates via input_stable.txt:
        createMake(make_template, POP_STRUC.POPULATION(Ind_No).ANGLES);
        % For debugging purposes:
        %disp(POP_STRUC.POPULATION(Ind_No).COORDINATES);
    end
    
    [nothing, nothing] = unix(['cat input.make0_' num2str(step) ' > input.make0']);
    [nothing, nothing] = unix(['cat input.key_'   num2str(step) ' > input.key']);
    
    
    POP_STRUC.POPULATION(Ind_No).numIons    = 'N/A';
    POP_STRUC.POPULATION(Ind_No).Vol        = 'N/A';
    POP_STRUC.POPULATION(Ind_No).symg       = 'N/A';
    POP_STRUC.POPULATION(Ind_No).S_order    = 'N/A';
    POP_STRUC.POPULATION(Ind_No).order      = 'N/A';
    POP_STRUC.POPULATION(Ind_No).struc_entr = 'N/A';
    POP_STRUC.POPULATION(Ind_No).K_POINTS   = 'N/A';
end

