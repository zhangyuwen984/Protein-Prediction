function Write_CP2K(Ind_No)


% USPEX Version 8.3.2
% Change: created

global POP_STRUC
global ORG_STRUC

numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;
atomType    = ORG_STRUC.atomType;

try
    [nothing, nothing] = unix(['cat cp2k_options_' num2str(Step) ' > cp2k.inp ']);
catch
    disp(['cp2k_options is not present for step ' num2str(Step)]);
    quit
end

Lattice = latConverter(LATTICE);
Lattice = Lattice';
Lattice(4:6) = Lattice(4:6)*180/pi;

[nothing, nothing] = unix('echo  \&SUBSYS  > subsys.uspex');
[nothing, nothing] = unix('echo     \&CELL >> subsys.uspex');
[nothing, nothing] = unix(['echo       ABC [angstrom] ' num2str(Lattice(1:3), 11) ' >> subsys.uspex']);
[nothing, nothing] = unix(['echo       ALPHA_BETA_GAMMA [deg] ' num2str(Lattice(4:6), 11) ' >> subsys.uspex']);

if     ORG_STRUC.dimension == 3
       [nothing, nothing] = unix(['echo       PERIODIC XYZ >> subsys.uspex ']);
elseif ORG_STRUC.dimension == 2
       [nothing, nothing] = unix(['echo       PERIODIC XY >> subsys.uspex ']);
elseif ORG_STRUC.dimension == 0
       [nothing, nothing] = unix(['echo       PERIODIC NONE >> subsys.uspex ']);
end

[nothing, nothing] = unix('echo     \&END >> subsys.uspex');
[nothing, nothing] = unix('echo     \&COORD >> subsys.uspex');
[nothing, nothing] = unix('echo        SCALED >> subsys.uspex');

coordLoop = 1;
for i = 1 : length(numIons)
 for j = 1 : numIons(i)
  [nothing, nothing] = unix(['echo ' megaDoof(atomType(i)) ' ' num2str(COORDINATES(coordLoop,:), 11) ' >> subsys.uspex']);
  coordLoop = coordLoop + 1;
 end
end

%%  print the Pressure 
[nothing, nothing] = unix([' echo EXTERNAL_PRESSURE \[kbar\] ' num2str(ORG_STRUC.ExternalPressure*10) ' > pressure.uspex' ]);
%%

[nothing, nothing] = unix('echo     \&END >> subsys.uspex');
[nothing, nothing] = unix('echo  \&END SUBSYS >> subsys.uspex');
%[nothing, nothing] = unix('echo  \&END FORCE_EVAL >> subsys.uspex');
