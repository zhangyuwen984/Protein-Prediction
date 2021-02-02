function Heredity_311(Ind_No)

global POOL_STRUC
global ORG_STRUC
global OFF_STRUC

info_parents = struct('parent',{}, 'fracFrac', {},'dimension', {},'offset', {}, 'enthalpy', {});
Molecule = struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CREATING Offspring with heredity %%%%%%%%%%%%%%%%
searching = 1;
count = 1;
while searching
    count = count + 1;
    if count > 50
       %disp('failed to do Heredity in 50 attempts, switch to Random');
       USPEXmessage(507,'',0);
       Random_311(Ind_No);
       break;
    end

    %choose one set of parents
    same = 1;
    % make sure you don't choose twice the same structure for this heredity
    while same
        par_one = find (ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
        par_two = find (ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
        num1 = POOL_STRUC.POPULATION(par_one(end)).numMols;
        num2 = POOL_STRUC.POPULATION(par_two(end)).numMols;
        if (par_one(end) ~= par_two(end)) & (sum(num1)>1) & (sum(num2)>1)
            same = 0;
            IND1 = par_one(end);
            IND2 = par_two(end);
        end
    end

    goodHeritage = 0;
    goodLattice = 0;
    goodComposition = 0;
    securityCheck = 0;
    while goodHeritage + goodLattice + goodComposition  ~=3

        securityCheck = securityCheck+1;
        [numMols, numBlocks, potentialOffspring, potentialLattice,fracFrac,dimension,offset, table]= heredity_molecule(IND1, IND2);
        
        goodComposition = CompositionCheck(numBlocks);
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
                      center = POOL_STRUC.POPULATION(IND1).MOLECULES(table(indMol,2)).MOLCENTER;
                      MOL    = POOL_STRUC.POPULATION(IND1).MOLECULES(table(indMol,2)).MOLCOORS;
                 else
                      center = POOL_STRUC.POPULATION(IND2).MOLECULES(table(indMol,2)).MOLCENTER;
                      MOL    = POOL_STRUC.POPULATION(IND2).MOLECULES(table(indMol,2)).MOLCOORS;
                 end                          

                 absMUT_COORD = potentialOffspring (indMol,:) * potentialLattice;
                 displacement = absMUT_COORD - center;

                 for inder=1:3
                     MOL(:,inder)=MOL(:,inder)+displacement(:,inder);
                 end
                     Molecule(indMol).MOLCOORS=MOL;

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

        if goodHeritage + goodLattice + goodComposition == 3
            OFF_STRUC.POPULATION(Ind_No).MOLECULES = Molecule;
            OFF_STRUC.POPULATION(Ind_No).LATTICE = potentialLattice;
            OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
            OFF_STRUC.POPULATION(Ind_No).numMols = numMols;
            OFF_STRUC.POPULATION(Ind_No).typesAList = typesAList;
            OFF_STRUC.POPULATION(Ind_No).MtypeLIST = MtypeLIST;
            OFF_STRUC.POPULATION(Ind_No).numBlocks = numBlocks;

            ID1 = POOL_STRUC.POPULATION(IND1).Number;
            ID2 = POOL_STRUC.POPULATION(IND2).Number;
            E1  = POOL_STRUC.POPULATION(IND1).enthalpy;
            E2  = POOL_STRUC.POPULATION(IND2).enthalpy;
            info_parents(1).parent = num2str([ID1 ID2]);
            info_parents.enthalpy = 0;
            fracFrac = [0 fracFrac 1];
            for i = 2:length(fracFrac)
               ratio=fracFrac(i)-fracFrac(i-1);
               if mod(i,2)==1
                  info_parents.enthalpy = info_parents.enthalpy+E1*ratio;
               else
                  info_parents.enthalpy = info_parents.enthalpy+E2*ratio;
               end
            end
            info_parents.enthalpy=info_parents.enthalpy;
            info_parents.fracFrac=fracFrac;
            info_parents.dimension=dimension;
            info_parents.offset=offset;
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

%%%%%%%%%%%% END CREATING Offspring with heredity %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
