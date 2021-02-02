function [candidate, newLattice, errorS] = symope_2D(nsym, numIons, volume, minD, thickness)

global ORG_STRUC

spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];

numIons0 = numIons;

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
fprintf(fp,'%4d   ! space group \n', nsym);
fprintf(fp,'%6.3f ! area of primitive unit cell\n', volume);
fprintf(fp,'1.0 1.0 90.0 ! vector of primitive unit cell\n');
fprintf(fp,'%6.3f ! thickness of layer\n', thickness);
fprintf(fp,'%4d   ! number of types of atoms\n', length(numIons1));

for i=1:length(numIons1)
    fprintf(fp,'%4d ', numIons1(i));
end

fprintf(fp, '! number of atoms of each type\n');

for ii = 1 : size(minD,1)
    for jj = 1 : size(minD,2)
        fprintf(fp, '%5.3f ', minD(ii,jj));
    end
end
fprintf(fp, '! minimum distance between atoms \n');
fprintf(fp, '0  ! symmetry operation \n');

if ORG_STRUC.fixRndSeed>0
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);



[a, b] = unix([spgBINDIR '/random_2d < rc.in > rc.out']);

candidate = zeros(sum(numIons),3);

if a == 0
    % read rc.out
    try
        handle = fopen('rc.out');
        f1 = ceil(length(numIons1)^2/10);
        for LL = 1 : 10+f1
            tmp = fgetl(handle); % repeats input file
        end
        lat_tmp = fgetl(handle); % lattice string - cell parameters:  36.84031  36.84031  36.84031  90.00000  90.00000  90.00000
        %  if isempty(findstr(lat_tmp, 'cell parameters'))
        %  lat_tmp = fgetl(handle); % lattice string - cell parameters:  36.84031  36.84031  36.84031  90.00000  90.00000  90.00000
        %  end
        error1 = findstr(lat_tmp, 'error');
        if isempty(error1)
            lat_1 = sscanf(lat_tmp,'%*s %*s %g %g %g %g %g %g');
            tmp = fgetl(handle); % repeats input file
            tmp = fgetl(handle); % repeats input file
            for LL = 1 : sum(numIons)
                atom_tmp = fgetl(handle); % lattice strings in Bohr units
                atom_1 = sscanf(atom_tmp,'%*s %*s %g %g %g');
                candidate(LL,1) = atom_1(1); candidate(LL,2) = atom_1(2); candidate(LL,3) = atom_1(3);
            end
            status = fclose(handle);
            lat_1(4:6) = lat_1(4:6)*pi/180;
            lat = latConverter(lat_1);
            new = 1;
        else
            status = fclose(handle);
            new = 0;
            lat = [1 0 0; 0 1 0; 0 0 1];
        end
    catch
        status = fclose(handle);
        new = 0;
        lat = [1 0 0; 0 1 0; 0 0 1];
    end
    
    if new
        errorS = 0;
    else
        errorS = 1;
    end
    
    newLattice = lat;
    
else
    status = ['Error while generating the crystal with symmetry group ' num2str(nsym) ' and number of atoms ' num2str(numIons)];
    errorStatus1 = a;
    errorStatus2 = b;
    newLattice = [1 0 0; 0 1 0; 0 0 1];
    errorS = 1;
end

if sum(numIons0) ~= size(candidate,1)
    errorS = 1;
end

