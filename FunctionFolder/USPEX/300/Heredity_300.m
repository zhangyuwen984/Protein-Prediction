function Heredity_300(Ind_No)

% Change: added multiparents heredity
% We will reselect parents, if USPEX cannot produce good Offspring within 100 attempts
% Last updated by Qiang Zhu (2013/10/16)

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

info_parents = struct('parent',{}, 'fracFrac', {},'dimension', {},'offset', {},'enthalpy', {});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CREATING Offspring with heredity %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;
searching = 1;
while searching
    count = count + 1;
    if count > 50 | length(POP_STRUC.POPULATION) == 1
        %disp('failed to do Heredity in 50 attempts, switch to Random');
        USPEXmessage(507,'',0);
        Random_300(Ind_No);
        break;
    end
    
    % choose one set of parents
    % We select randomly a dimension used for the spacial criteria of the heredity
    dimension = RandInt(1,1,[1,3]);
    
    if (ORG_STRUC.manyParents == 0) | (ORG_STRUC.manyParents > 1)
        same = 1;
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
        parents = zeros(100,1);
        for i = 1 : 50
            parents(2*i) = IND1;
            parents(2*i-1) = IND2;
        end
        
    else
        parents_rank = randperm(max(ORG_STRUC.tournament));
        parents = zeros(length(ORG_STRUC.tournament),1);
        parents_ind = zeros(length(POP_STRUC.ranking),1);
        j1 = 1;
        for i1 = 1:max(ORG_STRUC.tournament)
            par = find (ORG_STRUC.tournament > (parents_rank(i1)-1));
            if parents_ind(POP_STRUC.ranking(par(end))) == 0
                parents_ind(POP_STRUC.ranking(par(end))) = 1;
                parents(j1) = POP_STRUC.ranking(par(end));
                j1 = j1 + 1;
            end
        end
        % Take into account maxDistHeredity
        if (ORG_STRUC.doFing)&(POP_STRUC.generation > ORG_STRUC.doFing-1)
            for i1 = 2 : length(ORG_STRUC.tournament)
                badParent = 0;
                for j1 = 1 : i1-1
                    f1 = POP_STRUC.POPULATION(POP_STRUC.ranking(parents(i1))).FINGERPRINT;
                    f2 = POP_STRUC.POPULATION(POP_STRUC.ranking(parents(j1))).FINGERPRINT;
                    if cosineDistance(f1, f2, ORG_STRUC.weight) > ORG_STRUC.maxDistHeredity
                        badParent = 1;
                    end;
                end;
                if badParent
                    c1 = parents(i1);
                    for j1 = i1 : length(ORG_STRUC.tournament)-1
                        parents(j1) = parents(j1+1);
                    end;
                    parents(length(ORG_STRUC.tournament)) = c1;
                end;
            end;
        end;
    end
    
    % two parents give many slabs,
    % 2 - every slab is chosen independently form previous,
    % 3 - slabs are cut from cell with fixed offset
    
    %wait for suitable offspring (offspring fulfilling hard constraints)
    goodHeritage = 0;
    goodLattice = 0;
    securityCheck = 0;
    while goodHeritage + goodLattice ~=2
        
        securityCheck = securityCheck+1;
        offset=[];
        
        
        if ORG_STRUC.manyParents == 0
            if sum(POP_STRUC.POPULATION(IND1).numIons)==sum(POP_STRUC.POPULATION(IND2).numIons) % ORG_STRUC.minAt==ORG_STRUC.maxAt 
                [numIons, potentialOffspring, potentialLattice,fracFrac,dimension,offset]= heredity_final(IND1,IND2);
            else
                [numIons, numBlocks, potentialOffspring, potentialLattice,fracFrac,dimension,offset]= heredity_final_var(IND1,IND2);
            end
        else
            [numIons, potentialOffspring, potentialLattice,fracFrac,dimension,offset]= heredity_finalMP(parents, dimension);
        end
        
        % optimize the lattice
        if ~ORG_STRUC.constLattice
            coord = potentialOffspring*potentialLattice;
            [coord, potentialLattice] = optLattice(coord, potentialLattice);
            potentialOffspring = coord/potentialLattice;
        end
        
        goodHeritage = distanceCheck(potentialOffspring, potentialLattice, numIons, ORG_STRUC.minDistMatrice);
        goodHeritage = ( CompositionCheck(numIons/ORG_STRUC.numIons) &  goodHeritage );
        goodHeritage = ( goodHeritage & (sum(numIons)>0) ); % sometimes it will generate structures with no atoms, need to bo fixed here!
        goodLattice = latticeCheck(potentialLattice);
        
        if goodHeritage + goodLattice  == 2
            OFF_STRUC.POPULATION(Ind_No).COORDINATES = potentialOffspring;
            OFF_STRUC.POPULATION(Ind_No).LATTICE = potentialLattice;
            OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
            fracFrac = [0 fracFrac 1];
            enthalpy = 0;
            ID = [];
            for i = 2:length(fracFrac)
                ID= [ID POP_STRUC.POPULATION(parents(i-1)).Number];
                N_atom = sum(POP_STRUC.POPULATION(parents(i-1)).numIons);
                E = POP_STRUC.POPULATION(parents(i-1)).Enthalpies(end)/N_atom;
                ratio=fracFrac(i)-fracFrac(i-1);
                enthalpy = enthalpy+E*ratio;
            end
            info_parents(1).parent = num2str(ID);
            info_parents.enthalpy  = enthalpy(end);
            info_parents.enthalpy  = info_parents.enthalpy;
            info_parents.fracFrac  = fracFrac;
            info_parents.dimension = dimension;
            info_parents.offset    = offset;
            OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
            OFF_STRUC.POPULATION(Ind_No).howCome = ' Heredity ';
            disp(['Structure ' num2str(Ind_No) ' generated by heredity']);
            searching=0;
        end
        
        if securityCheck > 20
            % disp('Cannot produce good Offspring within 20 attempts')
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% END CREATING Offspring with heredity %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
