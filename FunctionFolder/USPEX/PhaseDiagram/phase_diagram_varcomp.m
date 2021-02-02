% Description: MatLab program to get phase diagram data from USPEX output (varcomp).
% Author     : Maxim Rakitin, maksim.rakitin@stonybrook.edu
% Date       : 2013-12-30
% Last update: 2014-02-17 - fixed problem when the list of ID's starts not
%                           from 1, but any other number.

% $Rev: 943 $
% $Author: maxim $
% $Date: 2015-03-08 22:58:06 +0400 (Sun, 08 Mar 2015) $

global ORG_STRUC
global USPEX_STRUC

%load('Current_ORG.mat');
%load('USPEX.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants:
eV_to_J   = 1.602176487E-19;
A3_to_m   = 1.0E-30;
Pa_to_GPa = 1.0E-9;
SI_factor = (eV_to_J * Pa_to_GPa/A3_to_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total number of structures/individuals:
tot_N = size(USPEX_STRUC.POPULATION, 2);
id_list = obtainDifferentStructures(tot_N);

% Get necessary data from USPEX_STRUC:
[comp_list, enth_list, volume_list] = parseIndividuals(id_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- This region is for input variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIons   = ORG_STRUC.numIons;
N_Block   = size(numIons,1);     %block_Type
N_Type    = size(numIons,2);     %atom_Type
N_Point   = length(id_list);
numBlocks = zeros(N_Point, N_Block);
Energy    = zeros(N_Point, 1);
for i=1:N_Point
    numBlocks(i,:) = USPEX_STRUC.POPULATION(id_list(i)).numBlocks;
    Energy(i)      = enth_list(i)/sum(numBlocks(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fitness, convex_hull] = phase_update_hull(N_Block, N_Point, numBlocks, Energy, numIons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execution and printing:
fid = fopen('convex_hull_pressures.txt', 'wb');

% 0 GPa:
X_console = sprintf('Reference pressure convex hull:');
disp(X_console);
fprintf(fid, strcat(X_console, '\n'));

format_header = '%16s %9s %9s %5s %5s %9s %7s';
format_data   = '[%14s] %9.4f %9.3f %5i %5i %9.3f %7.4f';
X_header_console = sprintf(format_header, 'Composition', 'Volume', 'Enthalpy', 'ID', 'Symm', 'Q_entropy', 'A_order');
X_header_file = strcat(X_header_console, '\n');

disp(X_header_console);
fprintf(fid, X_header_file);

for i=1:length(convex_hull(:,1))
    comp_str = num2str(sprintf('%s', num2str(uint8(convex_hull(i, 1:size(numIons, 1))))));

    id = convex_hull(i,end);
    X_console = sprintf(format_data,                           ...
                        comp_str,                              ...
                        volume_list(id),                       ...
                        enth_list(id),                         ...
                        id_list(id),                           ...
                        USPEX_STRUC.POPULATION(id).symg,       ...
                        USPEX_STRUC.POPULATION(id).struc_entr, ...
                   mean(USPEX_STRUC.POPULATION(id).order)      ...
                   );
    disp(X_console);
    fprintf(fid, strcat(X_console, '\n'));
end
disp(' ');
fprintf(fid, '\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -100 ... 100 GPa:
P0      = 0;     % Reference value
P_step  = 2;     % Step in GPa
counter = 0;
previous_ids      = [];
stable_POSCARS_ID = [];
for P=-100:P_step:100

    Energy    = zeros(N_Point, 1);
    for i=1:N_Point
        H0  = enth_list(i);
        V   = volume_list(i);
        H = H0 + V * (P - P0)/SI_factor;
        Energy(i) = H/sum(numBlocks(i,:));
    end
    
    [null, convex_hull] = phase_update_hull(N_Block, N_Point, numBlocks, Energy, numIons);
    convex_hull = sortrows(convex_hull, size(convex_hull, 2));
    
    current_ids = [];
    for i=1:size(convex_hull(:,end))
        current_ids = [current_ids, id_list(convex_hull(i,end))];
    end
    
    if isequal(sort(previous_ids), sort(current_ids)) == false
        X_console = sprintf('\n%i GPa pressure convex hull:', P);
        disp(X_console);
        fprintf(fid, strcat(X_console, '\n'));
        
        disp(X_header_console);
        fprintf(fid, X_header_file);
        for i=1:length(convex_hull(:,1))
            comp_str = num2str(sprintf('%s', num2str(uint8(convex_hull(i, 1:size(numIons, 1))))));
            id = convex_hull(i,end);
            X_console = sprintf(format_data,                           ...
                                comp_str,                              ...
                                volume_list(id),                       ...
                                enth_list(id),                         ...
                                id_list(id),                           ...
                                USPEX_STRUC.POPULATION(id).symg,       ...
                                USPEX_STRUC.POPULATION(id).struc_entr, ...
                           mean(USPEX_STRUC.POPULATION(id).order)      ...
                           );
            disp(X_console);
            fprintf(fid, strcat(X_console, '\n'));
        end
        counter = counter + 1;
        stable_POSCARS_ID = [stable_POSCARS_ID, current_ids];
    end
    previous_ids = current_ids;
end

X = sprintf('\n%i different convex hulls found.', counter);
disp(X);
fprintf(fid, X);

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\nThe following structures will be saved in stablePOSCARS_pressure file:'));
disp(num2str(stable_POSCARS_ID));
extractStableStructures(stable_POSCARS_ID);
disp(' ');
