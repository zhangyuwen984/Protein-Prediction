function Write_FHIaims(Ind_No)
global POP_STRUC
global ORG_STRUC
numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;

if exist('FHI_output')
   [nothing, nothing] = unix('mv FHI_output FHI_output_old');
end
if exist('relaxation_restart_file.FHIaims')
   [nothing, nothing] = unix('rm relaxation_restart_file.FHIaims');
end
if exist('geometry.in.next_step')
   %[nothing, nothing] = unix('mv geometry.in.next_step geometry_old');
end

try
   [nothing, nothing] = unix(['cat aims_control_' num2str(Step) ' > control.in']);
catch
   disp(['aims_control_ ' num2str(Step) 'NOT EXISTING']);
   quit
end

[Kpoints, Error] = Kgrid(LATTICE, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
if Error == 1  % This LATTICE is extremely wrong, let's skip it from now
    POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
    POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
end

if ORG_STRUC.dimension ~= 0
    [nothing, nothing] = unix(['echo k_grid ' num2str(Kpoints) ' >> control.in']);
end

Write_AIMS_Geometry(Ind_No);
