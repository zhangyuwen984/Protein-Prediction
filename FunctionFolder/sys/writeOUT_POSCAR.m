function writeOUT_POSCAR(Ind_No, doOrder)

global POP_STRUC
global ORG_STRUC
%need to rewrite this code in future
%Ind_No  = 0; for PSO best structure
%doOrder = 0; POSCAR for vasp calculation

if (Ind_No == 0)  % best structure in PSO regime
   lat     = POP_STRUC.PSO(POP_STRUC.bestPSOstruc).lattice;
   coord   = POP_STRUC.PSO(POP_STRUC.bestPSOstruc).coordinates;
   numIons = POP_STRUC.PSO(POP_STRUC.bestPSOstruc).numIons;
   order   = POP_STRUC.PSO(POP_STRUC.bestPSOstruc).order;
 bodyCount = 0;
else
   lat     = POP_STRUC.POPULATION(Ind_No).LATTICE;
   coord   = POP_STRUC.POPULATION(Ind_No).COORDINATES;
   numIons = POP_STRUC.POPULATION(Ind_No).numIons;
 bodyCount = POP_STRUC.POPULATION(Ind_No).Number;
   if doOrder
      order = POP_STRUC.POPULATION(Ind_No).order;
   end
end

atomType = ORG_STRUC.atomType;
SGtolerance = ORG_STRUC.SGtolerance;

[a,b]=unix(['cat /dev/null > POSCAR']);

Lattice_par = latConverter(lat);

if size(Lattice_par,1) == 6
   Lattice_par = Lattice_par';
end
Lattice_par(4:6) = Lattice_par(4:6)*180/pi;

if (ORG_STRUC.doSpaceGroup == 1)
   [nothing, current_path] = unix('pwd');
   current_path(end) = [];
   cd([ORG_STRUC.homePath '/CalcFoldTemp']);
   symg = anasym_stokes(lat, coord, numIons, atomType, SGtolerance);
   cd(current_path)
else
   symg = 0;
end

if Ind_No > 0  
   POP_STRUC.POPULATION(Ind_No).symg = symg;
end

if doOrder
  fp1 = fopen('POSCAR_order', 'w');
  fprintf(fp1, 'EA%-4d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f Sym.group: %4d\n', bodyCount, Lattice_par(1:6), symg);
  fprintf(fp1, '1.0000\n');
  for latticeLoop = 1 : 3
     fprintf(fp1, '%12.6f %12.6f %12.6f\n', lat(latticeLoop,:));
  end

     if ORG_STRUC.dimension==2
        Surface_numIons = POP_STRUC.POPULATION(Ind_No).Surface_numIons;
        Surface_order   = POP_STRUC.POPULATION(Ind_No).Surface_order;
        Surface_LATTICE = POP_STRUC.POPULATION(Ind_No).Surface_LATTICE;
        Surface_COOR = POP_STRUC.POPULATION(Ind_No).Surface_COORDINATES;
        Surface_COOR = Surface_COOR*Surface_LATTICE/lat;
        for i=1:length(Surface_numIons)
           fprintf(fp1, '%4d ', Surface_numIons(i));
        end
           fprintf(fp1, '\n');
           fprintf(fp1, 'Direct\n');
        for coordLoop = 1 : sum(Surface_numIons)
            fprintf(fp1, '%12.6f %12.6f %12.6f %10.4f\n', Surface_COOR(1,:), Surface_order(coordLoop));
        end
     else
       for i=1:length(numIons)
          fprintf(fp1, '%4d ', numIons(i));
       end
          fprintf(fp1, '\n');
          fprintf(fp1, 'Direct\n');
       for coordLoop = 1 : sum(numIons)
          fprintf(fp1, '%12.6f %12.6f %12.6f %10.4f\n', coord(coordLoop,:), order(coordLoop));
       end
     end
  fclose(fp1);
end

fp = fopen('POSCAR', 'w');
fprintf(fp, 'EA%-4d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f Sym.group: %4d\n', bodyCount, Lattice_par(1:6), symg);
fprintf(fp, '1.0000\n');

for latticeLoop = 1 : 3
   fprintf(fp, '%12.6f %12.6f %12.6f\n', lat(latticeLoop,:));
end

if (ORG_STRUC.varcomp == 1) & (doOrder==0)  %prepare for VASP POSCAR
   for i=1:length(numIons)
      if numIons(i)>0 
         fprintf(fp, '%4d ', numIons(i));
      end
   end
else
%   for i=1:length(numIons)
%      fprintf(fp,'%4s', megaDoof(ceil(ORG_STRUC.atomType(i))));
%   end
%      fprintf(fp, '\n');
   for i=1:length(numIons)
      fprintf(fp, '%4d ', numIons(i));
   end 
end
fprintf(fp, '\n');
fprintf(fp, 'Direct\n');

for coordLoop = 1 : sum(numIons)
    fprintf(fp, '%12.6f %12.6f %12.6f\n', coord(coordLoop,:));
end
fclose(fp);
