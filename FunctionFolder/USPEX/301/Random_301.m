function Random_301(Ind_No)

global ORG_STRUC
global OFF_STRUC
global POP_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

goodStructure = 0;
goodLattice = 0;


if exist('Seeds/compositions')
    %[a, b]=unix(['wc Seeds/compositions']); %Sometimes the composition file is broken
    %tmp = str2num(b(1:2));
    [a, b]=unix(['cat Seeds/compositions |wc -l']);
    tmp = str2num(b);
    if tmp(1)>0
        ORG_STRUC.firstGeneSplit=load('Seeds/compositions');
        ORG_STRUC.splitN = length(ORG_STRUC.firstGeneSplit);
    end
end


newSym = 1;
badLat   = 0;
badStruc = 0;
badSymmetry = 0;

minDistMatrice = ORG_STRUC.minDistMatrice;


tic
failedDist = 0;
while goodStructure + goodLattice  ~= 2
    failedTime = toc;
    if (failedDist > 10000) || (failedTime > 300)
        if minDistMatrice(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
            if failedTime > 300
                %disp('WARNING! Can not generate a structure after 5 minutes. The minimum distance threshold will be lowered by 0.1.');
                USPEXmessage(501,'',0);
            else
                %disp('WARNING! Can not generate a structure after 10000 tries. The minimum distance threshold will be lowered by 0.1.');
                USPEXmessage(502,'',0);
            end
            %disp('Please check your IonDistances parameter.');
            %disp(' ');
            minDistMatrice = 0.9*minDistMatrice;
            failedDist = 0;
            tic
        else
            disp('Could not generate a structure after 30000 tries or 15 minutes.');
            disp('Please check the input files. The calculation has to stop.');
            disp('Possible reasons:  unreasonably big IonDistances.');
            disp('but not too small for pseudopotential overlap errors to kill interatomic repulsion.')
            quit;
        end
    end
    
    errorS = 0;
    
    tmp = ceil(rand*ORG_STRUC.splitN);
    if tmp == 0
        tmp = 1;
    end
    numBlocks = ORG_STRUC.firstGeneSplit(tmp,:);
    numIons   = numBlocks*ORG_STRUC.numIons;
    candidate = rand(sum(numIons),3);
    %sym_coef = 1;
    
    if ORG_STRUC.constLattice    % lat1 = lattice
        lat1 = ORG_STRUC.lattice;
        lat = lat1;
    else                         % lat1 = lattice volume
        lat1 = 0;
        for it = 1 : length(ORG_STRUC.latVolume)
            lat1 = lat1 + numBlocks(it)*ORG_STRUC.latVolume(it);
        end
        new = 0;
        while ~new
            startLat = rand(6,1);
            startLat(4:6) = startLat(4:6)*(pi/2);
            check_startLat = latConverter(startLat);
            volLat = det(check_startLat);
            if volLat == 0
                continue;
            end
            % scale the lattice to the volume we asume it approximately to be
            ratio = lat1/volLat;
            startLat(1:3) = startLat(1:3)*(ratio)^(1/3);
            [dummy, lat2] = optLattice([0 0 0], latConverter(startLat)); % optimize lattice
            new = latticeCheck(lat2);
            if new
                lat = lat2;
            end
        end
    end
    
    sym_coef = ORG_STRUC.sym_coef;
    if (badSymmetry > 15) | ((badSymmetry > 5) & (sum(ORG_STRUC.splitInto)>3))
        badSymmetry = 0;
        newSym = 1;      % change the symmetry group if can't generate the crystal
    end
    badSymmetry = badSymmetry + 1;
    if newSym
        tmp = find(ORG_STRUC.nsym > 0);
        nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
        newSym = 0;
    end
    
    if nsym >= 1  % [MR]: it was > 1, changing it to >=1 since now we are able to generate P1 structures ccorrectly
        % OLD: if nsym == 1 - just keep the random candidate that we generated
        cd([ORG_STRUC.homePath '/CalcFoldTemp']);
        [candidate, lat, errorS] = symope_crystal(nsym, numIons, lat1, minDistMatrice, ORG_STRUC.sym_coef);
        
        if errorS == 0
            [resulted_sg_no, nothing] = determine_spacegroup(lat, candidate, ORG_STRUC.atomType, ...
                numIons, nsym, POP_STRUC.bodyCount + Ind_No, ...
                ORG_STRUC.SGtolerance, ORG_STRUC.homePath);
        end
        
        cd(ORG_STRUC.homePath)
    end
    
    if sum(ORG_STRUC.splitInto)>3 %splitcell
        [lat, errorS, candidate] = splitBigCell(latConverter(lat), ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand)),numIons, nsym);
    end
    
    if errorS == 0
        goodStructure = distanceCheck(candidate, lat, numIons, minDistMatrice*sym_coef);
        goodLattice = latticeCheck(lat);
        
        badLat = badLat + 1 - goodLattice;
        badStruc = badStruc + 1 - goodStructure;
        if badLat > 1000
            %disp('Lattice failed 1000 times during random structure generation.')
            %disp(' ')
            USPEXmessage(503,'',0);
            badLat = 0;
        end
        if badStruc > 1000
            %disp('Structure failed 1000 times during random structure generation.')
            %disp(' ')
            USPEXmessage(504,'',0);
            badStruc = 0;
        end
    else
        goodStructure = 0;
        goodLattice = 0;
    end
    
    if goodStructure + goodLattice == 2
        OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
        OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
        OFF_STRUC.POPULATION(Ind_No).Parents = [];
        OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
        OFF_STRUC.POPULATION(Ind_No).numBlocks = numBlocks;
        OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
        
        disp(['Structure ' num2str(Ind_No) ' built with the symmetry group  ' num2str(nsym) ' (' spaceGroups(nsym) ') composition ' num2str(numIons) ]);
    else
        failedDist = failedDist + 1;
    end
    
end

