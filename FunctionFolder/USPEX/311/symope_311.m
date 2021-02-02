function [candidate, newLattice, numSites, Operation, errorS, errorN, P, PB] = symope_311(nsym, numIons, lat, minD)


global ORG_STRUC

spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];

errorS = 0;
errorN = 0;
P  = [1 2 3];
PB = [1 2 3];
permutation = [1 2 3];   % no swaps
permutationBack = [1 2 3];
sGroup = spaceGroups(nsym); % space group's standard symbol
if (nsym<75) & (nsym>15) & (sGroup(1) == 'P')
    ListPerm = [1 2 3; 3 2 1; 1 3 2; 2 1 3; 3 1 2; 2 3 1];  %possible permutation;
    id = RandInt(1,1,6)+1;
    permutation = ListPerm(id,:); % all swaps allowed
    if id == 5
        permutationBack = ListPerm(6,:);  % exchange is not reversible (3 1 2) to (2 3 1)
    elseif id == 6
        permutationBack = ListPerm(5,:);  % exchange is not reversible (3 1 2) to (2 3 1)
    else
        permutationBack = permutation;
    end
end

if (size(lat,1) == 1) & (size(lat,2) == 1)        % random lattice, lat =  volume
    fixLat = 0;
    lat1 = zeros(1,6);
    lat1(1:3) = rand(1,3) + 0.5;
    x = 0;
    while (x < 0.3)
        lat1(4:6) = (rand(1,3)*120 + 30)*pi/180;% fortran works with degrees
        x = 1 - cos(lat1(4))^2 - cos(lat1(5))^2 - cos(lat1(6))^2 + 2*cos(lat1(4))*cos(lat1(5))*cos(lat1(6));
    end
    ratio = lat/det(latConverter(lat1));
    lat1(4:6) = lat1(4:6)*180/pi; % fortran works with degrees
    lat1(1:3) = lat1(1:3)*(ratio)^(1/3);
else                       % fixed lattice, lat = lattice
    fixLat=1;
    if (size(lat,1) == 3) & (size(lat,2) == 3)
        % swap axes
        %%% HAVE TO SWAP COLUMNS TOO!!!
        lat1 = latConverter(lat);
        lat2 = lat1;
        for axis = 1 : 3
            lat1(axis) = lat2(permutation(axis));
        end
    else
        % swap axes
        lat1 = lat;   % INPUT is 1*6
        for axis = 1 : 3
            lat1(axis) = lat(permutation(axis));
        end
    end
    lat1(4:6) = lat1(4:6)*180/pi; % fortran works with degrees
    if permutation == [3 2 1]  %if it is monoclinic, we also swap alpha and beta.
        lat2=lat1;
        lat1(4)=lat2(6);
        lat1(6)=lat2(4);
    end
end

i = 1;
while i <= length(numIons)
    if numIons(i) == 0
        numIons(i) = [];
        minD(i,:) = [];
        minD(:,i) = [];
        i = i - 1;
    end
    i = i + 1;
end
fp = fopen('rc.in', 'w');
fprintf(fp,'%4d ! space group \n', nsym);
fprintf(fp,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ! lattice of primitive unit cell\n', lat1(1:6));
fprintf(fp,'%4d ! number of types of atoms\n', length(numIons));

for i=1:length(numIons)
    fprintf(fp,'%4d ', numIons(i));
end
fprintf(fp, '! number of atoms of each type\n');

for ii = 1 : size(minD,1)
    for jj = 1 : size(minD,2)
        fprintf(fp, '%5.3f ', minD(ii,jj));
    end
end
fprintf(fp, '! minimum distance between atoms \n');
if ORG_STRUC.fixRndSeed==1
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);


[a, b] = unix([spgBINDIR '/random_cell_mol < rc.in > rc.out']);
[a, b]=unix(['grep NUMBER rc.out |wc -l']);
sites = str2num(b);

candidate = zeros(sum(numIons),3);
numSites = zeros(sites,length(numIons));
if a == 0
    %         FIX THE PROGRAM WHEN THERE IS AN ERROR IN THE OUTPUT FILE : instead of lattice - not possible to place given number of atoms into unit cell of given spac e-group symmetry
    % read rc.out
    handle = fopen('rc.out');
    % fixing possible multiline minDist (only 10 numbers per output line!)
    f1 = ceil(length(numIons)^2/10);
    for LL = 1 : 6 + f1
        tmp = fgetl(handle); % repeats input file
    end
    sumatom = 0;
    for i = 1 : sites
        tmp = fgetl(handle); % repeats input file
        tmp = fgetl(handle);
        if ~isempty(str2num(tmp))
            numSites(i,:) = str2num(tmp);
            if i>1
                numSites(i,1) = sum(numSites(i,:)) - sumatom;
            end
            sumatom = sumatom +numSites(i,1);
        end
    end
    
    [a, b]=unix(['grep operation rc.out |wc -l']);
    numo = str2num(b);
    Operation = zeros(numo,3);
    for i=1:numo
        a = fgetl(handle);
        Operation(i,1:3) = sscanf(a, '%*s %g %g %g', [1,3]);
    end
    lat_tmp = fgetl(handle); % lattice string - cell parameters:  36.84031  36.84031  36.84031  90.00000  90.00000  90.00000
    error1 = findstr(lat_tmp, 'error');
    if isempty(error1)
        lat_1 = sscanf(lat_tmp,'%*s %*s %g %g %g %g %g %g');
        lat_2 = zeros(1,6);
        lat_2(1) = lat_1(1); lat_2(2) = lat_1(2); lat_2(3) = lat_1(3);
        lat_2(4) = lat_1(4)*pi/180; lat_2(5) = lat_1(5)*pi/180; lat_2(6) = lat_1(6)*pi/180;
        tmp = fgetl(handle); % repeats input file
        tmp = fgetl(handle); % repeats input file
        for LL = 1 : sum(numIons)
            atom_tmp = fgetl(handle); % lattice strings in Bohr units
            atom_1 = sscanf(atom_tmp,'%*s %*s %g %g %g');
            candidate(LL,1) = atom_1(1); candidate(LL,2) = atom_1(2); candidate(LL,3) = atom_1(3);
        end
        status = fclose(handle);
        % try to fix the lattice first, if needed
        lat = latConverter(lat_2);
        %   abs_cand = candidate*lat_2;
        %   [abs_cand,lat] = optLattice(abs_cand,lat_2); % optimize lattice
        %   candidate = abs_cand/lat;
        
        new = latticeCheck(lat);
    else
        if findstr(lat_tmp, 'not possible')
            errorN = 1;
        end
        
        status = fclose(handle);
        new = 0;
        lat = [1 0 0; 0 1 0; 0 0 1];
        
    end
    
    newLattice = lat;
    if ~new
        errorS = 1;
    end
else
    status = ['Error while generating the crystal with symmetry group ' num2str(nsym) ' and number of atoms ' num2str(numIons)];
    newLattice = [1 0 0; 0 1 0; 0 0 1];
    errorS = 1;
end

% swap axes (and coordinates) back
if errorS == 0
    newLattice1 = latConverter(newLattice);
    newLattice2 = newLattice1;
    candidate1 = candidate;
    for axis = 1 : 3
        newLattice2(axis) = newLattice1(permutationBack(axis));
        candidate(:,axis) = candidate1(:,permutationBack(axis));
    end
    if permutation == [3 2 1]  %if it is monoclinic, we also swap alpha and beta.
        newLattice2(4)=newLattice1(5);
        newLattice2(5)=newLattice1(4);
    end
    newLattice = latConverter(newLattice2);
    P = permutation;
    PB = permutationBack;
end

