function extendedConvexHull_301(convex_hull, fit_max)

global ORG_STRUC
global USPEX_STRUC
% This function creates a list of structures that are closer than fit_max to convex hull
% read enthalpies and compositions
atomType = ORG_STRUC.atomType;
numIons  = ORG_STRUC.numIons;
popSize  = max(ORG_STRUC.populationSize, 20);  %in the case of a very small size

N        = length(USPEX_STRUC.POPULATION);
generation  = length(USPEX_STRUC.GENERATION);
[Nb, Ntype] = size(numIons);
enth    = zeros(N,1);
fitness = zeros(N,1);
comp    = zeros(N,Ntype);
comp1   = zeros(N,Nb);
for i=1:N
    comp(i,:)  = USPEX_STRUC.POPULATION(i).numIons;    %count by each atomType
    comp1(i,:) = USPEX_STRUC.POPULATION(i).numBlocks;  %count by each blocktype
    enth(i)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(comp(i,:));
end
% here we build the new fitness, using the convex hull
fitness = final_convex_hull(convex_hull, enth, comp1, numIons);
for i=1:N
    USPEX_STRUC.POPULATION(i).Fitness = fitness(i);
end

[nothing, ranking] = sort(fitness);
[accepted, newInd, Na] = RemoveStruc(comp, fitness, fit_max, ranking, USPEX_STRUC, generation);
%--------------------------------------------------------------------------
% Energy_base: the energy of each ending member in the system
Energy_base = zeros(size(comp1,2),1);
for i = 1: size(comp1,2)
    Energy_base(i)=convex_hull(i,end-1)*sum(numIons(i,:)); %ev/Block
end

%Initialization of some variables
x = zeros(sum(accepted(:)),size(numIons,2)-1);   
y = zeros(sum(accepted(:)),1);
x_= zeros(sum(newInd(:)),size(numIons,2)-1);
y_= zeros(sum(newInd(:)),1);
a1 = 0;
a1_new= 0;
count = 0; %count the structures of POOL

% output all accepted structures
FileName = 'extended_convex_hull';
fp1 = fopen(FileName, 'w');
[nothing, nothing] = unix(['cat /dev/null > extended_convex_hull_POSCARS']);  %Start to get new POSCAR
fprintf(fp1,'# It contains the all the information in extendedConvexHull.pdf\n');
fprintf(fp1,'# X axis: Composition \n');
fprintf(fp1,'# Y axis: Formation energy relative to the substance\n');
fprintf(fp1,'# Fitness: its distance to the convex hull (eV/block)\n');
fprintf(fp1,'  ID   Compositions    Enthalpies     Volumes     Fitness   SYMM    X        Y     \n');
if size(comp,2) == 2  %only two type of atoms, we build the convex_hull of AB-B
    fprintf(fp1,'                       (eV/atom)    (A^3/atom)   (eV/block)              (eV/atom) \n');
else
    fprintf(fp1,'                       (eV/atom)    (A^3/atom)   (eV/block)              (eV/block) \n');
end

for i = 1 : Na
    j = ranking(i);
    if accepted(j)>0
        a1 = a1 + 1;
        [x(a1,:), y(a1)] = Get_XY(numIons, comp(j,:), Energy_base, enth(j));
        if newInd(j)==1
           a1_new=a1_new+1;
           x_(a1_new,:)=x(a1,:);
           y_(a1_new)=y(a1);
        end
        symg      = USPEX_STRUC.POPULATION(j).symg;
        num       = USPEX_STRUC.POPULATION(j).numIons;
        lat       = USPEX_STRUC.POPULATION(j).LATTICE;
        coords    = USPEX_STRUC.POPULATION(j).COORDINATES;
        order     = USPEX_STRUC.POPULATION(j).order;
        numBlocks = USPEX_STRUC.POPULATION(j).numBlocks;
        volume    = USPEX_STRUC.POPULATION(j).Vol/sum(num);       
        composition = sprintf('%3d',num);
        shift=[4, 2, 1]; %so far we only consider 6 component
        if size(composition,2)<11
            composition=[composition,blanks(shift(length(num)))];
        end
        if size(numIons,2)==2
           fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d  %6.3f  %7.4f\n', ...
            j, composition, enth(j), volume, fitness(j), symg, x(a1), y(a1));
        elseif size(numIons,2)==3
           fprintf(fp1,'%4d  [%11s]   %9.4f    %9.4f   %9.4f   %4d  %6.3f %6.3f  %7.4f\n', ...
            j, composition, enth(j), volume, fitness(j), symg, x(a1,:), y(a1));
        end
        %POSCARS
        Write_POSCAR(atomType, j, symg, num, lat, coords);
        [nothing, nothing] = unix([' cat POSCAR      >> ' FileName '_POSCARS']);    
        %Build a POOL_STRUC
        if count <  popSize
    	   count = Update_POOL_STRUC(count, j, enth, fitness, USPEX_STRUC.POPULATION(j));        
        end 
    end
end
fclose(fp1);

% make convex hull figures here
makeFigure_convexhull_301(generation, convex_hull, comp1, [x'; y'], [(x_)';( y_)'], Energy_base);

