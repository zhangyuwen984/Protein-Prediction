function makeFigures_M400(pickUpNCount, Nsteps, ConstLat)
%function makeFigures_M400()

% $Rev: 1120 $
% $Author: mrakitin $
% $Date: 2015-08-20 21:51:42 +0400 (Thu, 20 Aug 2015) $

% All the data are loaded from USPEX_STRUC, and OUTPUT.txt and origin files.

global USPEX_STRUC

%pickUpNCount = 0;
%Nsteps       = 1;
%ConstLat     = 0;
%load USPEX.mat;

N         = length(USPEX_STRUC.POPULATION);
enth      = zeros(N,Nsteps);
vol       = zeros(N,1);
fit       = zeros(N,1);
ID        = zeros(N,1);

for i=1:N
    numIons    = USPEX_STRUC.POPULATION(i).numIons;
    fit(i)     = USPEX_STRUC.POPULATION(i).Fitness;
    enth(i,:)  = USPEX_STRUC.POPULATION(i).Enthalpies/sum(numIons);
    vol(i) = det(USPEX_STRUC.POPULATION(i).LATTICE)/sum(numIons);
    ID(i)  = i;
end


% creates a set of figures in tif format
% Nsteps = number of optimization steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    [nothing, nothing] = unix(['rm Fitness_vs_N.tif']);
    x = []; y = [];
    x = ID;
    y = fit;
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');

    % Define ranges +/-std from mean value:
    mean_y   = mean(y);
    std_y    = std(y);
    coef_std = 1;
    y_min    = mean_y - coef_std * std_y;
    y_max    = mean_y + coef_std * std_y;
    axis([min(x), max(x), y_min, y_max]);

    xlabel('Structure number');
    ylabel('Fitness');
    print(h,'-dtiff','-r120','Fitness_vs_N.tif');
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 (Echild - Eparent(s)) as function of N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nothing, nothing] = unix(['rm Variation-Operators.tif']);
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    % 0 - rotation, 1 - heredity
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
            if ~isempty(findstr(tmp, 'Rotation'))
                c(i) = 0;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Heredity'))
                c(i) = 1;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'SecSwitch'))
                c(i) = 2;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'ShiftBorder'))
                c(i) = 3;
                s(i) = 30;
            end
            CE = str2num(tmp(17:26));
            PE = str2num(tmp(27:36));
            y(i) = CE-PE;
        end
    end
    status = fclose(handle);
    subplot(2,1,1);
    
    % Define ranges +/-std from mean value:
    mean_y   = mean(y);
    std_y    = std(y);
    coef_std = 1;
    y_min    = mean_y - coef_std * std_y;
    y_max    = mean_y + coef_std * std_y;

    scatter(x,y,s,c,'filled','MarkerEdgeColor','k');
    xlabel('Structure index');
    ylabel('E_{Child} - E_{parent}');
    axis([min(x), max(x), y_min, y_max]);
    hcb = colorbar('YTick',[0:3],'YTickLabel',{'Rotation', 'Heredity', 'SecSwitch', 'ShiftBorder'});
    set(hcb,'YTickMode','manual')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fraction versus N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a, b] = unix(['grep frac OUTPUT.txt |cut -c60-65']);
    if ~isempty(b)
        x=[]; y=[];
        frac = str2num(b);
        N = size(frac,1)/5;   % Was 7
        x = 1:1:N;
        y = zeros(5,N);       % Was 7
        str = {'Heredity', 'Random', 'Rotation', 'SecSwitch', 'ShiftBorder'};
        for i = 1:N
            for j = 1:5       % Was 7
                y(j,i) = frac(5*(i-1)+j);     % Was 7
            end
        end
        subplot(2,1,2); 
        plot(x,y(1,:),x,y(2,:),x,y(3,:),x,y(4,:),x,y(5,:));
        legend(str);
        xlabel('Generation number');
        ylabel('Fraction of each variation');
        if N>1
            axis([1 N 0 1]);
        end
    end

    print(h,'-dtiff','-r120','Variation-Operators.tif');
catch
end
