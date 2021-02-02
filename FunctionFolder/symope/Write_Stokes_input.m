function Write_Stokes_input(numIons, minD, nsym, lat1, fixRndSeed, sym_coef);
numIons1 = numIons;
i = 1;
while i <= length(numIons1)
    if numIons1(i) == 0
        numIons1(i) = [];
        minD(i,:) = [];
        minD(:,i) = [];
        i = i - 1;
    end
    i = i + 1;
end

fp = fopen('rc.in', 'w');
fprintf(fp,'%4d ! space group \n', nsym);
fprintf(fp,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ! lattice of primitive unit cell\n', lat1(1:6));
fprintf(fp,'%4d ! number of types of atoms\n', length(numIons1));

for i=1:length(numIons1)
    fprintf(fp,'%4d ', numIons1(i));
end
fprintf(fp, '! number of atoms of each type\n');

minD = minD*sym_coef;
for ii = 1 : size(minD,1)
    for jj = 1 : size(minD,2)
        fprintf(fp, '%5.3f ', minD(ii,jj));
    end
end
fprintf(fp, '! minimum distance between atoms \n');
sym_coef = 1;
fprintf(fp, '%2d coefficient between minDist and symmetrization distance \n', sym_coef);
if fixRndSeed>0
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);
