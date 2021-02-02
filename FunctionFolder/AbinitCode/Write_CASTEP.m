function Write_CASTEP(Ind_No)
% created by Zamaan Raza

global POP_STRUC
global ORG_STRUC

atomType    = ORG_STRUC.atomType;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;

try
  [nothing, nothing] = unix(['cat param_' num2str(Step) ' > cstp.param']);
catch
  disp(['param file is not present for step' num2str(Step)]);
  quit
end

try
  [nothing, nothing] = unix(['cat cell_' num2str(Step) ' > cstp.cell']);
catch
  disp(['cell file is not present for step' num2str(Step)]);
  quit
end


%LATTICE
fp = fopen('cstp.cell', 'a+');
fprintf(fp, '%%BLOCK LATTICE_CART\n');
fprintf(fp, '%6.3f %6.3f %6.3f\n', LATTICE(1,:));
fprintf(fp, '%6.3f %6.3f %6.3f\n', LATTICE(2,:));
fprintf(fp, '%6.3f %6.3f %6.3f\n', LATTICE(3,:));
fprintf(fp, '%%ENDBLOCK LATTICE_CART\n');

%COORDINATES
fprintf(fp, '%%BLOCK POSITIONS_FRAC\n');
coordLoop = 1;
for i = 1 : length(numIons)
  for j = 1 : numIons(i)
    fprintf(fp, '%4s  %8.4f %8.4f %8.4f\n', megaDoof(atomType(i)), COORDINATES(coordLoop,:));
    coordLoop = coordLoop + 1;
  end
end
fprintf(fp, '%%ENDBLOCK POSITIONS_FRAC\n');

%%%%Kgrid
[Kpoints, Error] = Kgrid(LATTICE, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
if Error == 1  % This LATTICE is extremely wrong, let's skip it from now
    POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
    POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
end

fprintf(fp, 'KPOINT_MP_GRID %2d %2d %2d\n', Kpoints(1,:));
fclose(fp);
