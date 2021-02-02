% Description: MatLab program to prepare phase diagram data from USPEX output.
% Author     : Maxim Rakitin, maksim.rakitin@stonybrook.edu
% Date       : 2013-10-24
% Updated    : 2013-12-30

% $Rev: 685 $
% $Author: maxim $
% $Date: 2014-10-24 10:30:20 +0400 (Fri, 24 Oct 2014) $

global USPEX_STRUC

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

pressures_list = [];
for i=1:length(id_list)
    pressures_list = [pressures_list; 0];
end

init_matrix_IVEP = [id_list,     ...
                    volume_list, ...
                    enth_list,   ...
                    pressures_list];

init_matrix_IVEP = sortrows(init_matrix_IVEP, 2);   % Sorted by volume.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate phase transforms for negative pressures (volume increases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_matrix_IVEP = init_matrix_IVEP;
neg_phase_transforms = [];
while true
    matrix_IVEP = current_matrix_IVEP;

    if matrix_IVEP(1, 4) == 0
        pos = 1;
        matrix_IVEP = sortrows(matrix_IVEP, 3);    % Sorted by enthalpy (asc).
    else
        pos = size(matrix_IVEP, 1);
        matrix_IVEP = sortrows(matrix_IVEP, 4);    % Sorted by pressure (asc).
    end

    id        = matrix_IVEP(pos, 1);
    vol       = matrix_IVEP(pos, 2);
    enthalpy  = matrix_IVEP(pos, 3);
    pressure  = matrix_IVEP(pos, 4);

    neg_phase_transforms = [neg_phase_transforms; id, vol, enthalpy, pressure];

    current_matrix_IVEP = [];
    for i=1:size(matrix_IVEP, 1)
        if (matrix_IVEP(i, 3) > enthalpy) && (matrix_IVEP(i, 2) > vol)
            current_pressure = -(matrix_IVEP(i, 3) - enthalpy)* SI_factor / (matrix_IVEP(i, 2) - vol);
            current_matrix_IVEP = [current_matrix_IVEP; matrix_IVEP(i, 1), ...
														matrix_IVEP(i, 2), ...
														matrix_IVEP(i, 3), ...
														current_pressure  ...
								   ];
        end
    end

    if size(current_matrix_IVEP, 1) == 0
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate phase transforms for positive pressures (volume decreases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_matrix_IVEP = init_matrix_IVEP;
pos_phase_transforms = [];
while true
    matrix_IVEP = current_matrix_IVEP;

    if matrix_IVEP(1, 4) == 0
        pos = 1;
        matrix_IVEP = sortrows(matrix_IVEP, 3);    % Sorted by enthalpy (asc).
    else
        pos = 1;
        matrix_IVEP = sortrows(matrix_IVEP, 4);    % Sorted by pressure (asc).
    end
    
    id        = matrix_IVEP(pos, 1);
    vol       = matrix_IVEP(pos, 2);
    enthalpy  = matrix_IVEP(pos, 3);
    pressure  = matrix_IVEP(pos, 4);
    
    pos_phase_transforms = [pos_phase_transforms; id, vol, enthalpy, pressure];
    
    current_matrix_IVEP = [];
    for i=1:size(matrix_IVEP, 1)
        if (matrix_IVEP(i, 3) > enthalpy) && (matrix_IVEP(i, 2) < vol)
            current_pressure = -(matrix_IVEP(i, 3) - enthalpy)* SI_factor / (matrix_IVEP(i, 2) - vol);
            current_matrix_IVEP = [current_matrix_IVEP; matrix_IVEP(i, 1), ...
														matrix_IVEP(i, 2), ...
														matrix_IVEP(i, 3), ...
														current_pressure  ...
								   ];
        end
    end
    
    if size(current_matrix_IVEP, 1) == 0
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos_phase_transforms = sortrows(pos_phase_transforms, 3);
all_values_by_volume = sortrows([pos_phase_transforms([2:end],:); neg_phase_transforms], 2); % Sort by volume

ids      = all_values_by_volume(:, 1);
volume   = all_values_by_volume(:, 2);
enthalpy = all_values_by_volume(:, 3);
pressure = all_values_by_volume(:, 4);

% Print figure with all values and found "pressure" lines:
figure('Name','Pressures output', 'NumberTitle','off'); %, 'Visible','off');
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing

plot(...
     volume, enthalpy,...
     'Color','blue',...
     'LineStyle',':',...
     'Marker','o',...
     'MarkerSize',8.0,...
     'MarkerEdgeColor',[0.5,0.5,0.5],...
     'MarkerFaceColor','c'...
    );
title('Pressures output');
xlabel('Volume, A^3/atom');
ylabel('Enthalpy, eV/atom');

labels = num2str(ids);  %' # ID labels of points
text(volume, enthalpy, labels,...
     'VerticalAlignment','middle', ...
     'HorizontalAlignment', 'right',...
     'BackgroundColor',[.7 .9 .7],...
     'Margin',0.5,...
     'FontSize',12,...
     'Rotation',90 ...
    )

% Add all volume-enthaply values:
hold on
plot(init_matrix_IVEP(:,2), init_matrix_IVEP(:,3),...
     'LineStyle','None',...
     'Marker','o',...
     'MarkerFaceColor','green',...
     'MarkerEdgeColor','black',...
     'MarkerSize',4.0...
    );
hold off
print -dpng Pressures_output.png;


% Save text data to a file:
fid = fopen('Pressures_output.txt', 'wb');
stable_POSCARS_ID = [];

X_file = sprintf('%12s %12s %9s %5s %3s %5s %6s\n', 'Volume', 'Enthalpy', 'Pressure', 'ID', 'SYMM', 'Q_entropy', 'A_order');  % Organize header
fprintf(fid, X_file);
X_console = sprintf('%12s %12s %9s %5s %3s %5s %6s', 'Volume', 'Enthalpy', 'Pressure', 'ID', 'SYMM', 'Q_entropy', 'A_order');  % Organize header
disp(X_console);
for i=1:size(all_values_by_volume, 1)
    X_file = sprintf('%12f %12f %9.1f %5i %4i %9.3f %7.4f\n', volume(i), ...
															  enthalpy(i), ...
															  pressure(i), ...
															  ids(i), ...
															  USPEX_STRUC.POPULATION(ids(i)).symg, ...
															  USPEX_STRUC.POPULATION(ids(i)).struc_entr, ...
														 mean(USPEX_STRUC.POPULATION(ids(i)).order) ...
					);
    fprintf(fid, X_file);
    X_console = sprintf('%12f %12f %9.1f %5i %4i %9.3f %7.4f', volume(i), ...
					   										   enthalpy(i), ...
															   pressure(i), ...
															   ids(i), ...
															   USPEX_STRUC.POPULATION(ids(i)).symg, ...
															   USPEX_STRUC.POPULATION(ids(i)).struc_entr, ...
														  mean(USPEX_STRUC.POPULATION(ids(i)).order) ...
					);
    disp(X_console);
    stable_POSCARS_ID = [stable_POSCARS_ID; ids(i)];
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\nThe following structures will be saved in stablePOSCARS_pressure file:'));
disp(num2str(stable_POSCARS_ID'));

extractStableStructures(stable_POSCARS_ID);

disp(' ');
