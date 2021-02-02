function [goodPop, POP, nsym] = Random_Init_310(goodPop, numMols0, nsym)


global ORG_STRUC


%Here we check the position of molecule center is reasonable or not
badSymmetry = 0;
newSym = 1;

CenterminDist = ORG_STRUC.CenterminDistMatrice;
numMols = numMols0; %ORG_STRUC.numMols;
running = 1;
while running
    if badSymmetry > 5
        badSymmetry = 0;
        newSym = 1;   % change the symmetry group if can't generate the crystal
    end
    badSymmetry = badSymmetry + 1;
    if newSym
        tmp = find(ORG_STRUC.nsym > 0);
        nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
        newSym = 0;
    end
    if ORG_STRUC.constLattice    % lat1 = lattice
        lat1 = ORG_STRUC.lattice;
    else                         % lat1 = lattice volume
        lat1 = ORG_STRUC.latVolume*sum(numMols)/sum(ORG_STRUC.numMols);
    end
    if (nsym > 1)  % if nsym == 1 - just keep the random candidate that we generated
        splitInto=ORG_STRUC.splitInto;
        if sum(splitInto)>3
            numMols=numMols0/prod(splitInto);
            for i=1:3
                if size(lat1)==[1 1]
                    lat1=lat1/splitInto(i);
                elseif size(lat1)==[3 3]
                    lat1(i,:)=lat1(i,:)/splitInto(i);
                else
                    lat1(i)=lat1(i)/splitInto(i);
                end
            end
        end
        cd(['CalcFoldTemp']);
        [candidate, lat, numSites, Operation, errorS, errorN, P, PB] = symope_310(nsym, numMols, lat1, CenterminDist);
        cd(ORG_STRUC.homePath);
        if errorN == 1
            ORG_STRUC.nsym(nsym)=0;
            newSym = 1;
        elseif errorS == 0
            POP.LATTICE = lat;
            running = 0;
            if distanceCheck(candidate, lat, numMols, CenterminDist-0.2)
                goodbad = 1;
            else
                % disp('Perhaps a bug of symmetry Code');
                % candidate
                goodbad = 0;
            end
        end
    end
end

if goodbad
    [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
    for item=1:20
        Molecules = GetOrientation(candidate, lat, numSites(:,1), Operation, MtypeLIST, nsym, P, PB);
        goodBad = newMolCheck(Molecules,lat, MtypeLIST, ORG_STRUC.minDistMatrice);
        if goodBad
            if sum(splitInto)>3
                [Molecules, numMols, lat] = SuperMol(Molecules, numMols, lat, splitInto);
                [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
            end
            
            [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
            POP.MOLECULES  = Molecules;
            POP.numMols    = numMols;
            POP.MtypeLIST  = MtypeLIST;
            POP.typesAList = typesAList;
            POP.numIons    = numIons;
            POP.LATTICE    = lat;
            POP.howCome    = '  Random  ';
            

            goodPop = goodPop+1;
            disp(['Molecular Crystal '  num2str(goodPop) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym)  ')']);
            
            break
        end
    end
end
