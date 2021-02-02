function Heredity_310(Ind_No)

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

info_parents = struct('parent',{},'fracFrac', {},'dimension', {},'offset', {}, 'enthalpy', {});
Molecule = struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{},'order',{},'Operation',{},'P',{},'PB',{});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CREATING Offspring with heredity %%%%%%%%%%%%%%%%
searching = 1;
count = 1;
while searching
    count = count + 1;
    if count > 50 | length(POP_STRUC.POPULATION) == 1
       %disp('failed to do Heredity in 50 attempts, switch to Random');
       USPEXmessage(507,'',0);
       Random_310(Ind_No);
       break;
    end

    %choose one set of parents
    same = 1;
    % make sure you don't choose twice the same structure for this heredity
    tryTime = 0;
    while same
        par_one = find (ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
        par_two = find (ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
        tryTime = tryTime + 1;
        if tryTime > 100   % To prevent the dead loop when too small popsize
           par_one = 1;
           par_two = 2;
        end
        if (par_one(end) ~= par_two(end))
            same = 0;
            IND1 = POP_STRUC.ranking(par_one(end));
            IND2 = POP_STRUC.ranking(par_two(end));
        end
    end

    goodHeritage = 0;
    goodLattice = 0;
    securityCheck = 0;
    while goodHeritage + goodLattice ~=2
        securityCheck = securityCheck+1;
        if norm(POP_STRUC.POPULATION(IND1).numMols-POP_STRUC.POPULATION(IND2).numMols)< 1e-3
            [numMols, potentialOffspring, potentialLattice,fracFrac,dimension,offset, table]= heredity_molecule(IND1,IND2);
            numBlocks = numMols/ORG_STRUC.numMols;
            goodComposition = 1;
        else
            [numMols, numBlocks, potentialOffspring, potentialLattice,fracFrac,dimension,offset, table]= heredity_molecule_var(IND1, IND2);
            goodComposition = CompositionCheck(numBlocks);
        end
        
        if (sum(numBlocks(:)==0)>0)
            goodLattice=0;
        else
            if ORG_STRUC.constLattice
                goodLattice=1;
            else
               goodLattice = latticeCheck(potentialLattice);
               coord = potentialOffspring*potentialLattice;
               [coord, potentialLattice] = optLattice(coord, potentialLattice);
               potentialOffspring = coord/potentialLattice;
            end
        end

        if goodLattice && goodComposition
            Molecule=[];
            for indMol= 1: sum(numMols)
                if table(indMol,1)==1      
                     center = POP_STRUC.POPULATION(IND1).MOLECULES(table(indMol,2)).MOLCENTER;
                     MOL    = POP_STRUC.POPULATION(IND1).MOLECULES(table(indMol,2)).MOLCOORS;
                else
                     center = POP_STRUC.POPULATION(IND2).MOLECULES(table(indMol,2)).MOLCENTER;
                     MOL    = POP_STRUC.POPULATION(IND2).MOLECULES(table(indMol,2)).MOLCOORS;
                end                          

                absMUT_COORD = potentialOffspring (indMol,:) * potentialLattice;
                displacement = absMUT_COORD - center;

                Molecule(indMol).MOLCOORS=bsxfun(@plus, MOL, displacement);
                
                if(size(MOL,1)==1)
                     Molecule(indMol).ZMATRIX=MOL;
                else
                     Molecule(indMol).ZMATRIX=real(NEW_coord2Zmatrix(MOL,ORG_STRUC.STDMOL(table(indMol,3)).format));
                end
            end
               [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);              
               goodHeritage = newMolCheck(Molecule, potentialLattice,MtypeLIST, ORG_STRUC.minDistMatrice-0.2);
        else
            goodHeritage=0;
        end

        if goodHeritage + goodLattice ==2
            OFF_STRUC.POPULATION(Ind_No).MOLECULES = Molecule;
            OFF_STRUC.POPULATION(Ind_No).LATTICE = potentialLattice;
            OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
            OFF_STRUC.POPULATION(Ind_No).numMols = numMols;
            OFF_STRUC.POPULATION(Ind_No).typesAList = typesAList;
            OFF_STRUC.POPULATION(Ind_No).MtypeLIST = MtypeLIST;
            parents = [IND1 IND2];
            fracFrac = [0 fracFrac 1];
            enthalpy = 0;
            ID = [];
            for i = 2:length(fracFrac)
                ID= [ID POP_STRUC.POPULATION(parents(i-1)).Number];
                E = POP_STRUC.POPULATION(parents(i-1)).Enthalpies(end);
                ratio=fracFrac(i)-fracFrac(i-1);
                enthalpy = enthalpy+E*ratio;
            end
            info_parents(1).parent = num2str(ID);
            info_parents.enthalpy  = enthalpy(end);
            info_parents.enthalpy  = info_parents.enthalpy/sum(numIons);
            info_parents.fracFrac  = fracFrac;
            info_parents.dimension = dimension;
            info_parents.offset    = offset;
            OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
            OFF_STRUC.POPULATION(Ind_No).howCome = ' Heredity ';
            disp(['Structure ' num2str(Ind_No) '  generated by heredity']);
            searching=0;
        end

        if securityCheck>20
            %  loop =loop-1;
            %we won't wait for a good offspring forever, will we?
            break
        end
    end
end

