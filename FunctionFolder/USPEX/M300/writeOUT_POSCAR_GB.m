function writeOUT_POSCAR_GB(Ind_No)

% USPEX Version 7.3.5
global POP_STRUC

lat = POP_STRUC.POPULATION(Ind_No).LATTICE;
coord = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
bodyCount = POP_STRUC.bodyCount;
Lattice_par = latConverter(lat);
Lattice_par(4:6) = Lattice_par(4:6)*180/pi;

symg='no SG';

fp = fopen('POSCAR', 'w');
fprintf(fp, 'EA%-4d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f Sym.group: %4s\n', bodyCount, Lattice_par(1:6), symg);
fprintf(fp, '1.0000\n');
  for latticeLoop = 1 : 3
     fprintf(fp, '%12.6f %12.6f %12.6f\n', lat(latticeLoop,:));
  end

  for i=1:length(numIons)
     fprintf(fp, '%4d ', numIons(i));
  end
     fprintf(fp, '\n');
     fprintf(fp, 'Selective dynamics\n');
     fprintf(fp, 'Direct\n');

 for i=1:size(coord,1)
     if POP_STRUC.POPULATION(Ind_No).chanAList(i)==1
        fprintf(fp, '%12.6f %12.6f %12.6f  T  T  T\n', coord(i,:));
     else
        fprintf(fp, '%12.6f %12.6f %12.6f  F  F  F\n', coord(i,:));
     end
 end

