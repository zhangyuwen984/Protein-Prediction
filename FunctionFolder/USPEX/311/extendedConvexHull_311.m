function extendedConvexHull_311(convex_hull, fit_max, min_dist)

% consider [1 0], [1,1], plot [5 1] [2 1], if we consider eV/atom
% 0-1/6-1/3-1/2
% 0-1/5-1/2-1 ratio will change
% This function creates a list of structures that are closer than fit_max to convex hull
% read enthalpies and compositions
% Octave cannot plot figure, let's skip it
% Last updated by Qiang Zhu (2013/12/27)

global POOL_STRUC
global ORG_STRUC
global USPEX_STRUC

% This function creates a list of structures that are closer than fit_max to convex hull
% read enthalpies and compositions
atomType = ORG_STRUC.atomType;
stopCrit = ORG_STRUC.stopCrit;
popSize  = max(ORG_STRUC.populationSize, 20);  %in the case of a very small size
bestHM   = max(ORG_STRUC.keepBestHM, 4);       %in the case of a very small size
numIons  = ORG_STRUC.numMols;
weight   = ORG_STRUC.weight;


N        = length(USPEX_STRUC.POPULATION);
generation  = length(USPEX_STRUC.GENERATION);
[Nb, Ntype] = size(numIons);
enth    = zeros(N,1);
fitness = zeros(N,1);
comp    = zeros(N,Ntype);
comp1   = zeros(N,Nb);
for i=1:N
    comp(i,:)  = USPEX_STRUC.POPULATION(i).numMols;
    comp1(i,:) = USPEX_STRUC.POPULATION(i).numBlocks;
    enth(i)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(comp(i,:));
end

% here we build the new fitness, using the convex hull
fitness = final_convex_hull(Nb, N, convex_hull, enth, comp1, numIons);
for i=1:N
    USPEX_STRUC.POPULATION(i).Fitness = fitness(i);
end
[nothing, ranking] = sort(fitness);

% now remove all structures above threshold from consideration
accepted = zeros(1, N);
newInd   = zeros(1, N);
Na = 0; % number of accepted structures
for i = 1 : N
    if fitness(ranking(i)) <= fit_max
        Na = Na + 1;
        if fitness(ranking(i)) >= -0.05 %remove structures with wrong energy
            accepted(ranking(i)) = 1;
            if USPEX_STRUC.POPULATION(ranking(i)).gen == generation
                newInd(ranking(i)) = 1;
            end
        end
    else
        break;
    end
end


% Remove all duplicates:
%--------------------------------------------------------------------------
% First step: fine filtering:
for i = 2 : Na
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j = 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
                    if SameStructure(ranking(i), ranking(j), USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                        break;
                    end
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
% Second step: rough filtering using A_order:
for i = 2 : Na
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j = 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if sameComposition(comp(ranking(i),:),comp(ranking(j),:))
                    if SameStructure_order(ranking(i), ranking(j), USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                        accepted(ranking(i)) = 0;
                        break;
                    end
                end
            end
        end
    else
        accepted(ranking(i)) = 0;
    end
end
%--------------------------------------------------------------------------


% output all accepted structures
% CHANGED in 9.2.0: plot based on blocks, not atomic types!
% for this we changed comp into comp1 in the code, where it was relevant

% make figures, if there are 2 types of blocks (atoms - before 9.2.0)
if size(comp1,2) == 2
    A0=convex_hull(1,end-1); %ev/atom
    B0=convex_hull(2,end-1); %eV/atom
end

x = zeros(1, sum(accepted(:)));
y = zeros(1, sum(accepted(:)));

x_= zeros(1, sum(newInd(:)));
y_= zeros(1, sum(newInd(:)));

a1 = 0;
a1_new= 0;
count = 0;
fp1 = fopen('extended_convex_hull', 'w');
[nothing, nothing] = unix(['cat /dev/null > extended_convex_hull_POSCARS']);  %Start to get new POSCAR
fprintf(fp1,'# It contains the all the information in extendedConvexHull.pdf\n');
fprintf(fp1,'# X axis: X = y/(x+y) \n');
fprintf(fp1,'# Y axis: Y = (E(AxBy)-x*E(A)-y*E(B))/(x+y)\n');
fprintf(fp1,'# Fitness: its distance to the convex hull\n');
fprintf(fp1,'# Note that sometimes the unit of Y is eV/block\n');
fprintf(fp1,'  ID   Compositions    Enthalpies     Volumes     Fitness   SYMM    X        Y     \n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)   (eV/block)              (eV/block)\n');

for i = 1 : Na
    j = ranking(i);
    if accepted(j)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make figures, if there are 2 types of atoms
        if size(comp1,2) == 2
            a1 = a1 + 1;
            x(a1) = comp1(j,2)/(comp1(j,1) + comp1(j,2));
            y(a1) = (enth(j)*sum(comp(j,:)) - A0*sum(numIons(1,:))*comp1(j,1) ...
                - B0*sum(numIons(2,:))*comp1(j,2)) / sum(comp1(j,:)); % eV/block
            if newInd(j)==1
                a1_new=a1_new+1;
                x_(a1_new)=x(a1);
                y_(a1_new)=y(a1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        symg      = USPEX_STRUC.POPULATION(j).symg;
        num       = USPEX_STRUC.POPULATION(j).numIons;
        lat       = USPEX_STRUC.POPULATION(j).LATTICE;
        coords    = USPEX_STRUC.POPULATION(j).COORDINATES;
        order     = USPEX_STRUC.POPULATION(j).order;
        MOLECULES = USPEX_STRUC.POPULATION(j).MOLECULES;
        numBlocks = USPEX_STRUC.POPULATION(j).numBlocks;
        numMols   = USPEX_STRUC.POPULATION(j).numMols;
        volume    = USPEX_STRUC.POPULATION(j).Vol/sum(num);
        
        composition = sprintf('%3d',num);
        shift=[4, 2, 1]; %so far we only consider 6 component
        if size(composition,2)<11
            composition=[composition,blanks(shift(length(num)))];
        end
        
        if size(comp1,2) == 2
        	fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d  %6.3f  %6.3f\n', ...
            	j, composition, enth(j), volume, fitness(j), symg, x(a1), y(a1));
        else
        	fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d\n', ...
            	j, composition, enth(j), volume, fitness(j), symg);
        end
        
        % POSCARS
        Write_POSCAR(atomType, j, symg, num, lat, coords);
        [nothing, nothing] = unix([' cat POSCAR      >> extended_convex_hull_POSCARS']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %New! We build a POOL_STRUC
        if count <  popSize
            ratio = numMols/sum(numMols);
            store = 0;
            newcomposition=1;
            for ind=1:size(POOL_STRUC.Composition_ratio,1)
                if sum(abs(ratio-POOL_STRUC.Composition_ratio(ind,:))) < 0.01  % old composition
                    if POOL_STRUC.Composition_surviving(ind)<stopCrit %we don't allow
                        E_ranking=POOL_STRUC.Composition_ranking(ind) + 1;  %Composition_ranking, how many structures
                        if E_ranking == 1
                            enthalpy1 = enth(j);
                            POOL_STRUC.Composition_fitness(ind) = fitness(j);
                            if enthalpy1 < POOL_STRUC.Composition_Bestenthalpy(ind)
                                POOL_STRUC.Composition_Bestenthalpy(ind) = enthalpy1;   %update bestenthalpy
                                POOL_STRUC.Composition_surviving(ind) = 1;
                            else
                                POOL_STRUC.Composition_surviving(ind) = POOL_STRUC.Composition_surviving(ind)+1;
                            end
                        end
                        if E_ranking < bestHM
                            store = 1;
                            POOL_STRUC.Composition_ranking(ind) = E_ranking;
                        end
                    else
                        POOL_STRUC.Composition_ranking(ind) = -1 ; %means we do not want this composition
                    end
                    newcomposition = 0;
                    break;
                end
            end
            
            if newcomposition==1;  %new composition
                POOL_STRUC.Composition_ratio(end+1,:)= ratio;
                POOL_STRUC.Composition_ranking(end+1) = 1;
                POOL_STRUC.Composition_surviving(end+1)= 1;
                POOL_STRUC.Composition_fitness(end+1)= fitness(j);
                POOL_STRUC.Composition_Bestenthalpy(end+1)= enth(j);
                store = 1;
            end
            if store == 1  % we select only a few structures for each composition
                count     = count + 1;
                POOL_STRUC.POPULATION(count).LATTICE = lat;
                POOL_STRUC.POPULATION(count).COORDINATES = coords;
                POOL_STRUC.POPULATION(count).numIons = num;
                POOL_STRUC.POPULATION(count).numMols = numMols;
                POOL_STRUC.POPULATION(count).numBlocks = numBlocks;
                POOL_STRUC.POPULATION(count).order = order;
                POOL_STRUC.POPULATION(count).MOLECULES = MOLECULES;
                POOL_STRUC.POPULATION(count).enthalpy = enth(j);
                POOL_STRUC.POPULATION(count).Number = j;
                
            end
        end %count reaches maximum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
fclose(fp1);

% make convex hull figures here

makeFigure_convexhull_311(generation, convex_hull, comp1, [x; y], [x_; y_]);
% try
%     if size(numIons,1) == 2
%         % add point to convex hull array, for plotting
%         inConvexHull = zeros(2, size(convex_hull,1));
%         for i = 1:size(convex_hull,1)
%             N_block = convex_hull(i,1:2);
%             N_atom = N_block*numIons;
%             E = convex_hull(i,3);
%             inConvexHull(1,i) = N_block(2)/sum(N_block); % y/(x+y) by block
%             inConvexHull(2,i) = (E*sum(N_atom) - A0*sum(numIons(1,:))*N_block(1) - B0*sum(numIons(2,:))*N_block(2)) / sum(N_block); % eV/block
%         end
%         
%         [nothing, chRanking] = sort(inConvexHull(1,:)); % we want line to go from left to right :)
%         h1 = figure;
%         set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
%         scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
%         xlabel('Composition');
%         ylabel('Enthalpy of formation (eV/block)');
%         line(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
%         hold on
%         box on;
%         scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)),'MarkerEdgeColor','g','MarkerFaceColor','k');
%         scatter(x_(1:a1_new), y_(1:a1_new), 'MarkerEdgeColor','k','MarkerFaceColor','r');
%         %   print(h1,'-dtiff','-r300','extendedConvexHull.tif');
%         print(h1,'-dpdf' ,        'extendedConvexHull.pdf');
%         print(h1,'-dpdf' ,       ['generation' num2str(generation)  '/extendedConvexHull.pdf']);
%         hold off
%     end
% catch
% end


function [fitness] = final_convex_hull(N_T, N_P, convex_hull, enthalpies, blocks, numIons)

%N_T = number of different blocks
%N_P = number of structures
fitness = zeros(1, N_P);

% Calculating fitness = distance from the point to convex hull
% To check the compounds for decomposition, we use binary string where 1 means that structure is in the basis for decomposition
% those elments build matrix C with their atom coefficients, X decomposable if CX = A has positive vector as a solution and det(C) <> 0

for it = 1 : N_P
    
    A = blocks(it,:);
    A = A'; % make it a column
    compos = zeros(1, size(convex_hull,1));
    for i = 1 : N_T
        compos(i) = 1;
        C(i,:) = convex_hull(i, 1:N_T);
    end
    fitn = enthalpies(it);
    fitness(it) = 1000000;
    C = C';
    
    while 1 % convex_hull cycle
        
        if det(C) ~= 0
            %X = C\A; 
            X = pinv(C)*A;
            if sum(X<-0.001)==0 % all X component must be positive
                f = 0;
                k = 1;
                for i = 1 : size(convex_hull,1)
                    if compos(i)
                        f = f + X(k)*convex_hull(i, N_T + 1)*sum(convex_hull(i, 1:N_T)*numIons);
                        k = k + 1;
                    end
                end
                f = f/sum(blocks(it,:)*numIons);
                if (f - fitn) < fitness(it)
                    fitness(it) = f - fitn;
                end
            end
        end
        
        flag = 1;
        for i = size(convex_hull,1) - N_T + 1 : size(convex_hull,1)
            if compos(i) == 0
                flag = 0;
            end
        end
        if flag
            break;
        end
        % update compos
        k = size(convex_hull,1) - 1;
        s = 1;
        while ~((compos(k) == 1) & (compos(k+1) == 0))
            s = s + compos(k+1);
            k = k - 1;
            if k < 1
                break; % just for impossible case if previous cycle fails, for example length(convex_hull) < N_T
            end
        end
        compos(k) = 0;
        compos(k+1:k+s) = 1;
        compos(k+s+1 : size(convex_hull,1)) = 0;
        k = 1;
        for i = 1 : size(convex_hull,1)
            if compos(i)
                C(k,:) = convex_hull(i, 1:N_T);
                k = k + 1;
            end
        end
        C = C';
    end
    
end


fitness = -1*fitness;
% very important and often neglected step;
% make sure precision is e-4 since it's compared with fitness of that precision!
fitness = round(fitness*10000)/10000;
