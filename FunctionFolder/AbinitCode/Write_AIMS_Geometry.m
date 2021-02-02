function Write_AIMS_Structure(Ind_No)

global POP_STRUC
global ORG_STRUC

numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;

fp = fopen('geometry.in','w+');

if ORG_STRUC.dimension ~= 0

   fprintf(fp, 'lattice_vector   %12.6f %12.6f %12.6f\n',  LATTICE(1,:));
   fprintf(fp, 'lattice_vector   %12.6f %12.6f %12.6f\n',  LATTICE(2,:));
   fprintf(fp, 'lattice_vector   %12.6f %12.6f %12.6f\n',  LATTICE(3,:));
   if ORG_STRUC.dimension == 2
      fprintf(fp, 'constrain_relaxation .true.\n');
   end
end


direct_coord = COORDINATES*LATTICE;
coordLoop = 1;
for i = 1 : length(numIons)
    for j = 1 : numIons(i)
        fprintf(fp, 'atom   %12.6f %12.6f %12.6f  %4s\n', direct_coord(coordLoop,:), megaDoof(ORG_STRUC.atomType(i)));
        if ORG_STRUC.dimension == 2
           if direct_coord(coordLoop,3) < (ORG_STRUC.bulk_lat(3,3) - ORG_STRUC.thicknessB)
              fprintf(fp, 'constrain_relaxation .true.\n');
           end
        end
        coordLoop = coordLoop + 1;
    end
end

fclose(fp);

