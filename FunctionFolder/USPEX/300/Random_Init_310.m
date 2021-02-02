function [ newCoords, lat] = Random_Init_310(Ind_No, numMols)

% implemented - USPEX Version 9.3.7
global ORG_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newSym = 1;
badSymmetry = 0;
goodBad = 0;
Molecules = struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{});

while ~goodBad
    [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
    if ORG_STRUC.constLattice    % lat1 = lattice
        lat1 = ORG_STRUC.lattice;
    else                         % lat1 = lattice volume
        lat1 = ORG_STRUC.latVolume*sum(numIons)/sum(ORG_STRUC.numIons);
    end
    if badSymmetry > 15
        badSymmetry = 0;
        newSym = 1;      % change the symmetry group if can't generate the crystal
    end
    
    badSymmetry = badSymmetry + 1;
    
    if newSym
        tmp = find(ORG_STRUC.nsym > 0);
        nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
        newSym = 0;
    end
    CenterMat = ORG_STRUC.CenterminDistMatrice;
    
    if nsym > 1  % if nsym == 1 - just keep the random candidate that we generated
        cd([ORG_STRUC.homePath '/CalcFoldTemp']);
        [candidate, lat, numSites, Operation, errorS, errorN, P, PB] = symope_310(nsym, numMols, lat1, CenterMat);
        cd(ORG_STRUC.homePath)
        if errorN==1
            ORG_STRUC.nsym(nsym)=0;
            newSym = 1;
        elseif errorS ==0
            if distanceCheck(candidate, lat, numMols, CenterMat-0.2)
                goodBad = 1; %sometime the stokes code fails
            else
                goodBad = 0;
            end
        end
    end
    
    if goodBad
        for item=1:30
            Molecules = [];
            [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
            Molecules = GetOrientation(candidate, lat, numSites(:,1), Operation, MtypeLIST, nsym, P, PB);
            for innerInder= 1: length(MtypeLIST)
                format = ORG_STRUC.STDMOL(MtypeLIST(innerInder)).format;
                Molecules(innerInder).ZMATRIX = NEW_coord2Zmatrix(Molecules(innerInder).MOLCOORS,format);
            end
            goodBad = newMolCheck(Molecules,lat, MtypeLIST, ORG_STRUC.minDistMatrice);
            if goodBad
                if  ORG_STRUC.maxAt > ORG_STRUC.minAt
                    disp(['Structure ' num2str(Ind_No) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym)  '), fragment ' num2str(numMols) ]);
                else
                    disp(['Structure ' num2str(Ind_No) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym)  ')' ]);
                end
                cattedCoors = [];
                for j = 1: sum(numMols)
                    cattedCoors = cat(1,cattedCoors,Molecules(j).MOLCOORS);
                end
                saveded = cattedCoors/lat;
                newCoords = zeros(0,3);
                for m = 1: length(ORG_STRUC.atomType)
                    s = find(typesAList == ORG_STRUC.atomType(m));
                    newCoords = cat(1,newCoords, saveded(s,:));
                end
                break
            end
        end
    end
end

