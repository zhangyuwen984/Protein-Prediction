function extendedConvexHull_201(convex_hull, E_AB, bulk_stoi, fit_max, flag)

% This function creates a list of structures that are closer than fit_max to convex hull
% Octave cannot plot the figure, let's skip it
% a new way to define InConvexHull
% last updated by Qiang Zhu (2014/01/11)
global USPEX_STRUC

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

% read enthalpies and compositions
atomType= USPEX_STRUC.SYSTEM(1).atomType;
weight  = USPEX_STRUC.SYSTEM(1).Fp_weight;

N       = length(USPEX_STRUC.POPULATION);
Ntype   = 2; %length(USPEX_STRUC.POPULATION(1).numIons);
enth    = zeros(N,1);
fitness = zeros(N,1);
comp    = zeros(N,Ntype);
cell    = zeros(N,2);
for i=1:N
    comp(i,:) = USPEX_STRUC.POPULATION(i).Surface_numIons;
    cell(i,:) = USPEX_STRUC.POPULATION(i).cell;
    enth(i)   = USPEX_STRUC.POPULATION(i).Enthalpies(end);
end

if flag == 2
    m = bulk_stoi(1);
    n = bulk_stoi(2);
end
% transform enthalpies into formation enthalpy
% transform compositions into x-coordinates
X = zeros(N, 1);
E = zeros(N, 1);
np = size(convex_hull,1);
for i = 1 : N
    E_DFT = enth(i);
    N_A = comp(i,1);
    N_B = comp(i,2);
    N_cell = prod(cell(i,:));
    if flag == 1
        E(i) = (E_DFT - N_A*E_AB)/N_cell; %Here E_AB = E_A
        X(i) = N_B/N_cell;
    else
        E(i) = (E_DFT - N_B*E_AB/n)/N_cell;
        X(i) = (N_A-m/n*N_B)/N_cell;
    end
end
fitness = toconvexhull(convex_hull, E, X);
[nothing, ranking] = sort(fitness);

% now remove all structures above threshold from consideration
accepted = zeros(1, N);
Na = 0; % number of accepted structures
for i = 1 : N
    if fitness(ranking(i)) <= fit_max
        Na = Na + 1;
        if fitness(ranking(i))>= -0.05  %no point below the established hull
            accepted(ranking(i)) = 1;
        end
    else
        break;
    end
end

% remove all duplicates
for i = 2 : Na
    ID1 = ranking(i);
    if USPEX_STRUC.POPULATION(ID1).ToCount > 0
        f1 = USPEX_STRUC.POPULATION(ID1).FINGERPRINT;
        for j = 1 : i-1
            ID2 = ranking(j);
            if USPEX_STRUC.POPULATION(ID2).ToCount > 0
                if ( cell(ID1,:)==cell(ID2,:) ) & ( comp(ID1,:)==comp(ID2,:) )
                    f2 = USPEX_STRUC.POPULATION(ID2).FINGERPRINT;
                    dist = cosineDistance(f1, f2, weight);
                    diff = abs(enth(ID1)-enth(ID2))/sum(comp(ID1,:));
                    if (dist < 0.01) & (diff < 0.01)
                        USPEX_STRUC.POPULATION(ID1).ToCount = 0;
                        accepted(ID1) = 0;
                        break;
                    end
                end
            end
        end
    else
        accepted(ID1) = 0;
    end
end

% make figures, if there are 2 types of blocks (atoms - before 9.2.0)
x = zeros(1, sum(accepted(:)));
y = zeros(1, sum(accepted(:)));
y_old = zeros(1, sum(accepted(:)));
inConvexHull = zeros(2, np); % x, y for convex hull points
x01=convex_hull(1,1);
y01=convex_hull(1,2);
x02=convex_hull(np,1);
y02=convex_hull(np,2);
for i = 1:np
    inConvexHull(1,i) = convex_hull(i,1); % y/(x+y) by block
    inConvexHull(2,i) = convex_hull(i,2)-y02-(convex_hull(i,1)-x02)*(y02-y01)/(x02-x01); % E
end

fp1 = fopen('extended_convex_hull', 'w');
fprintf(fp1,'# It contains all the information in extendedConvexHull.pdf\n');
if flag == 1
    fprintf(fp1,'# X axis: Delta_N  = N_B/N_cell \n');
    fprintf(fp1,'# Y axis: Delta_E0 = (E_DFT - N_A*E_A/n)/N_cell \n');
else
    fprintf(fp1,'# X axis: Delta_N  = (N_A-m/n*N_B)/N_cell \n');
    fprintf(fp1,'# Y axis: Delta_E0 = (E_DFT - N_B*E_AB/n)/N_cell \n');
end
fprintf(fp1,'# Fitness: its distance to the convex hull\n');
fprintf(fp1,'# Note the data in extendedConvexHull.pdf has been transformed to a horizontal hull.\n');
fprintf(fp1,'# Refence: Zhu, Q, et al, PRB, 87, 195317, 2013\n');

fprintf(fp1,'  ID   Compositions   Cell   Delta_N    Delta_E0     Fitness \n');
fprintf(fp1,'                             (/cell)    (eV/cell)   (eV/cell)\n');
[nothing, nothing] = unix(['cat /dev/null > extended_convex_hull_POSCARS']);  %Start to get new POSCAR

a1 = 1;
iC = 1; % convex hull index
for i = 1 : N
    j = ranking(i);
    if accepted(j)
        num       = USPEX_STRUC.POPULATION(j).numIons;
        lat       = USPEX_STRUC.POPULATION(j).LATTICE;
        coords    = USPEX_STRUC.POPULATION(j).COORDINATES;
        
        composition = sprintf('%3d',num);
        shift=[4, 2, 1]; %so far we only consider 6 component
        if size(composition,2)<11
            composition=[composition,blanks(shift(length(num)))];
        end
        
        fprintf(fp1,'%4d  [%11s]  [%1d %1d] %8.4f %12.4f %10.4f\n', ...
            j, composition, cell(j,:), X(j), E(j), fitness(j));
        % POSCARS
        Write_POSCAR(atomType, j, 0, num, lat, coords);
        [nothing, nothing] = unix([' cat POSCAR      >> extended_convex_hull_POSCARS']);
        % enthalpy of formation per atom = [E(AxBy) - (x*E0(A)+y*E0(B)]/(x+y)
        x(a1) = X(j);
        y(a1) = E(j)-y02-(x(a1)-x02)*(y02-y01)/(x02-x01);
        a1 = a1 + 1;
    end
end

fclose(fp1);
[nothing, chRanking] = sort(inConvexHull(1,:)); % we want line to go from left to right :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%octave cannot plot figure, let's skip it
try
    h1 = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    
    xlabel('\delta(N)');
    ylabel('Enthalpy of formation (eV)');
    line(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
    hold on
    box on;
    
    if OctaveMode == 0
        scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)),'MarkerEdgeColor','g','MarkerFaceColor','k');
    else
        scatter(inConvexHull(1,chRanking(:)),inConvexHull(2,chRanking(:)));
    end
    
    print(h1,'-dpdf','extendedConvexHull.pdf');
    hold off
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  fit = toconvexhull(hull, E, X);

np = size(hull,1);
N = length(E);
for loop = 1:N
    X10 = min(hull(:,1));
    X20 = max(hull(:,1));
    done=0;
    for i = 1:np
        if abs(X(loop)-hull(i,1))<0.001  %stable point
            fit(loop) = E(loop)- hull(i,2);
            done = 1;
            break
        elseif X(loop) < hull(i,1)
            if hull(i,1) < X20 + 0.001
                X2 = i;
                X20 = hull(i,1);
            end
        elseif X(loop) > hull(i,1)
            if hull(i,1) > X10 - 0.001
                X1 = i;
                X10 = hull(i,1);
            end
        end
    end
    if ~done
        E0= hull(X1,2)+  (X(loop)-X10)*(hull(X1,2)-hull(X2,2))/(X10-X20);
        fit(loop)=E(loop)-E0;
    end
end
