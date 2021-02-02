function Write_MOPAC(Ind_No)
global POP_STRUC
global ORG_STRUC
if exist('input.mop');
    [nothing, nothing] = unix('rm -rf input.mop input.arc  input.den  input.res');
end
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
step = POP_STRUC.POPULATION(Ind_No).Step;
[nothing, nothing] = unix(['cp mop_' num2str(step) '  input.mop']);
direct_coord = POP_STRUC.POPULATION(Ind_No).COORDINATES*POP_STRUC.POPULATION(Ind_No).LATTICE;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
fp = fopen('input.mop','a+');
coordLoop = 1;
if ORG_STRUC.dimension == 0
        for i = 1 : length(numIons)
                for j = 1 : numIons(i)
                        fprintf(fp, '%4s   %12.6f %12.6f %12.6f \n' , megaDoof(ORG_STRUC.atomType(i)), direct_coord(coordLoop,:));
                        coordLoop = coordLoop + 1;
                end
        end
elseif ORG_STRUC.dimension == 3 
        for i = 1 : length(numIons)
                for j = 1 : numIons(i)
                        fprintf(fp, '%4s   %12.6f %12.6f %12.6f \n' , megaDoof(ORG_STRUC.atomType(i)), direct_coord(coordLoop,:));
                        coordLoop = coordLoop + 1;
                end
        end
	for  a = 1 : 3
	        fprintf(fp, '%4s   %12.6f %12.6f %12.6f\n', 'Tv' , lattice(a,:));
	end
end
fclose(fp);


