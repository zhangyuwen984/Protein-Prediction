function makeFigures(pickUpNCount, Nsteps, ConstLat)

%All the data are loaded from USPEX_STRUC
%Lastly updated by Qiang Zhu (2014/02/18)

global USPEX_STRUC

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

N         = length(USPEX_STRUC.POPULATION);
enth      = zeros(N,Nsteps);
vol       = zeros(N,1);
fit       = zeros(N,1);
ID        = zeros(N,1);

for i=1:N
    numIons   = USPEX_STRUC.POPULATION(i).numIons;
    fit(i)    = USPEX_STRUC.POPULATION(i).Fitness;
    enth(i,:) = USPEX_STRUC.POPULATION(i).Enthalpies/sum(numIons);
    vol(i) =det(USPEX_STRUC.POPULATION(i).LATTICE)/sum(numIons);
    ID(i) = i;
end


% creates a set of figures in pdf format
% Nsteps = number of optimization steps
% E(n), fitness(n), E(V), Echild(Eparents), En(En-1) where n is optimization step
% if ConstLat = 1, we don't output E(V).
% Octave is not able to plot this figure, so let's skip it
% Last updated by Qiang Zhu (2014/02/20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    [a,b]=unix(['rm Energy_vs_N.pdf']);
    x = []; y = [];
    x = ID;
    y = enth(:,Nsteps);
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    xlabel('Structure number');
    ylabel('Enthalpy');
    print(h,'-dpdf', 'Energy_vs_N.pdf');
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    [a,b]=unix(['rm Fitness_vs_N.pdf']);
    x = []; y = [];
    x = ID;
    y = fit;
    y_max = max(y(find(y<1000))); % this is to kick out all the structures with fitness=10000;
    y(find(y>=10000))=y_max;
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    xlabel('Structure number');
    ylabel('Fitness');
    print(h,'-dpdf','-r120','Fitness_vs_N.pdf');
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ConstLat~=1
    try
        [a,b]=unix(['rm Energy_vs_Volume.pdf']);
        x = []; y = [];
        x = enth(:,Nsteps);
        y = vol;
        h = figure;
        set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
        if OctaveMode == 0
            scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y);
        end
        xlabel('Volume');
        ylabel('Enthalpy');
        print(h,'-dpdf','Energy_vs_Volume.pdf');
    catch
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs En(En-1) - to compare the energy change between different optimisation steps
try
    e_complete = enth;
    sub = ceil((Nsteps-1)/2);
    [a,b]=unix(['rm E_series.pdf']);
    for g = 1 : Nsteps - 1
        x = []; y = [];
        x = e_complete(:,g);
        y = e_complete(:,g+1);
        subplot(sub,2,g);
        if OctaveMode == 0
            scatter(x,y,10,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y,10);
        end
        xlabel(['E' num2str(g)]);
        ylabel(['E' num2str(g+1)]);
    end
    print(h,'-dpdf', 'E_series.pdf');
catch
end

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1 (Echild - Eparent(s)) as function of N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a,b]=unix(['rm Variation-Operators.pdf']);
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    % 0 - latmutation,  1 - softmutation & coormutation, 2 - rotation,
    % 3 - permutation,  4 - heredity
    N_count = N-pickUpNCount;
    x = ID(N-N_count+1:N);
    y = zeros(N_count,1);
    c = zeros(N_count,1);
    s = ones(N_count,1);
    handle = fopen('origin');
    tmp = fgetl(handle);
    for i=1:N_count
        tmp = fgetl(handle);
        if tmp ~= -1
            if ~isempty(findstr(tmp, 'LatMutated'))
                c(i) = 0;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'CoorMutate'))  % merged in softmutation!!!
                c(i) = 1;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'softmutate'))
                c(i) = 1;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Rotated'))
                c(i) = 2;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Permutate'))
                c(i) = 3;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Heredity'))
                c(i) = 4;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Spinmutate'))
                c(i) = 5;
                s(i) = 30;
            end
            CE = str2num(tmp(17:26));
            PE = str2num(tmp(27:36));
            y(i) = CE-PE;
        end
    end
    status = fclose(handle);
    subplot(2,1,1);
    if OctaveMode == 0
        scatter(x,y,s,c,'filled','MarkerEdgeColor','k');
    else
        scatter(x,y,s,c);
    end
    xlabel('Structure index');
    ylabel('E_{Child} - E_{parent}');
    hcb = colorbar('YTick',[0:1:4],'YTickLabel',{'Latmut.','AtomMut.','Rotation','Permut.','Heredity','Spinmutate'});
    set(hcb,'YTickMode','manual')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fraction versus N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a, b] = unix(['grep frac OUTPUT.txt |cut -c60-65']);
    %numOperater = 8;
    if ~isempty(b)
        x=[]; y=[];
        frac = str2num(b);
        N = size(frac,1)/8;
        x = 1:1:N;
        y = zeros(8,N);
        str = {'Heredity', 'Random', 'Softmutation', 'Permutation', 'Latmutation', 'Rotation', 'Transmutation', 'Spinmutate'};
        for i = 1:N
            for j = 1:8
                y(j,i) = frac(8*(i-1)+j);
            end
        end
        subplot(2,1,2);
        plot(x,y(1,:),x,y(2,:),x,y(3,:),x,y(4,:),x,y(5,:),x,y(6,:),x,y(7,:),x,y(8,:));
        legend(str);
        xlabel('Generation number');
        ylabel('Fraction of each variation');
        if N>1
            axis([1 N 0 1]);
        end
    end
    
    print(h,'-dpdf', 'Variation-Operators.pdf');
catch
end
