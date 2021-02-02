function [candidate, lat] = Random_Init_300(Ind_No, numIons)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global POP_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodStructure = 0;
goodLattice = 0;

newSym = 1;
badSymmetry = 0;
failedDist = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
badStruc = 0;
tic

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
            disp('Remember they should be much smaller than the real interatomic distances,');
            disp('but not too small for pseudopotential overlap errors to kill interatomic repulsion.')
            quit;
        end
    end
    
    errorS = 0;
    candidate = rand(sum(numIons),3);
    sym_coef = 1;
    if ORG_STRUC.constLattice    % lat1 = lattice
        if ORG_STRUC.minAt ~= ORG_STRUC.maxAt
            error('Single block doesnot support constantLattice, Please check your INPUT.txt file')
            quit();
        end
        lat1 = ORG_STRUC.lattice;
        lat = lat1;
    else                         % lat1 = lattice volume
        lat1 = ORG_STRUC.latVolume*sum(numIons)/sum(ORG_STRUC.numIons);
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
    
    if ORG_STRUC.topologyRandom == 0
        if ORG_STRUC.nsymN(1,1) == 0  % H. Stokes code to create a crystal with given symmetry
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
            if sum(ORG_STRUC.splitInto) > 3 %splitcell
                [lat, errorS, candidate] = splitBigCell(latConverter(lat), ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand)),numIons, nsym);
            end
        end
    else
        lat_volume = ORG_STRUC.latVolume;
        
        str_latVolume = num2str(lat_volume);
        sGroup = 'NONE';%spaceGroups(nsym);
        str_numIons = num2str(numIons);
        str_Ind_No = num2str(Ind_No);
        coordinatingNumbers = ORG_STRUC.coordinatingNumbers;
        safesave ([ORG_STRUC.homePath '/coordinatingNumbers.mat'],  coordinatingNumbers);
        cd([ORG_STRUC.homePath '/FunctionFolder/topology']);
        result_top = python_uspex([ORG_STRUC.homePath '/FunctionFolder/USPEX/src/random_topology.py'],str_Ind_No,str_latVolume,sGroup,str_numIons,1);
        cd(ORG_STRUC.homePath);
        [nothing, nothing] = unix('rm coordinatingNumbers.mat');
        [nothing, nothing] = unix('rm coordinatingNumbers.mat.backup');
        %disp(num2str(result_top));
        topology_number = result_top(1);
        lat = transpose(reshape(result_top(2:10),3,3));
        candidate = transpose(reshape(result_top(11:end),3,[]));
    end
    
    if errorS == 0
        goodStructure = distanceCheck(candidate, lat, numIons, minDistMatrice*sym_coef);
        goodLattice = latticeCheck(lat);
        
        badStruc = badStruc + 1 - goodStructure;
        
        if badStruc > 1000
            %disp('Structure failed 1000 times during random structure generation.')
            %disp(' ')
            USPEXmessage(503,'',0);
            badStruc = 0;
        end
    else
        goodStructure = 0;
        goodLattice = 0;
    end
    
    if goodStructure + goodLattice == 2
        % Structures with wrong space group:
        if nsym >= 1
            if nsym ~= resulted_sg_no
                if ~isfield(ORG_STRUC, 'wrong_spacegroups')
                  ORG_STRUC.wrong_spacegroups = 1;
                else
                  ORG_STRUC.wrong_spacegroups = ORG_STRUC.wrong_spacegroups + 1;
                end
            end
        end
        
        if ORG_STRUC.maxAt > ORG_STRUC.minAt
            if ORG_STRUC.topologyRandom == 0
                disp(['Structure ' num2str(Ind_No) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym) ') , composition ' num2str(numIons)]);
            else
                disp(['Structure ' num2str(Ind_No) ' built with the topology ' num2str(topology_number) ' , composition ' num2str(numIons) ]);
            end
        else
            if ORG_STRUC.topologyRandom == 0
                disp(['Structure ' num2str(Ind_No) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym) ')']);
            else
                disp(['Structure ' num2str(Ind_No) ' built with the topology ' num2str(topology_number)  ]);
            end
        end
    else
        failedDist = failedDist + 1;
    end
    
end
