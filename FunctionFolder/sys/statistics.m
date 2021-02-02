function statistics(target_enthalpy)
% $Rev: 1006 $
% $Author: qian $
% $Date: 2015-04-25 15:36:13 +0400 (Sat, 25 Apr 2015) $

%% Function to walk in the specified path recursively:
function dirs = getAllDirs(path, dirs)
    listing = dir(path);

    for i = 1:size(listing,1)
        if listing(i).isdir & strcmp(listing(i).name, '..') == 0 & ...
                strcmp(listing(i).name, '.') == 0
            dirs{end+1} = [path '/' listing(i).name];
            dirs = getAllDirs([path '/' listing(i).name], dirs);
        end
    end
end

%% Whole list of directories:
dirs = getAllDirs(pwd, {});

%% Filter results* directories where results* is a basedir
results_dirs = {};
for i = 1:size(dirs,2)
    dirs{i};
    if strfind(dirs{i}, 'results') > 0
        [pathstr, name, ext] = fileparts(dirs{i});
        if strfind(name, 'results') > 0
            %disp(dirs{i});
            results_dirs{end+1} = dirs{i};
        end
    end
end

%% Check if USPEX.mat exists in the found results* folders to perform analysis:
USPEX_mat = {};
for i = 1:size(results_dirs,2)
    if exist([results_dirs{i} '/' 'USPEX.mat'], 'file')
        %disp([results_dirs{i} '/' 'USPEX.mat']);
        USPEX_mat{end+1} = [results_dirs{i} '/' 'USPEX.mat'];
    end
end


%% Analyze found USPEX.mat files:
number_of_mat = size(USPEX_mat,2);

disp(' ');
disp(['Number of files to be processed: ' num2str(number_of_mat)]);
disp(['Target enthalpy: ' num2str(target_enthalpy)]);
disp(' ');

if number_of_mat > 0
    found_generations = -1*ones(number_of_mat,1);
    found_populations = -1*ones(number_of_mat,1);
    found_enthalpies  = nan(number_of_mat,1);

    for i = 1:number_of_mat
        %disp(USPEX_mat{i})

        % Load USPEX_STRUC from found .mat file:
        global USPEX_STRUC
        load(USPEX_mat{i});

        % Get enthalpies from populations:
        pop_enth_full = nan(size(USPEX_STRUC.POPULATION(:)));
        pop_enth      = [];
        for j=1:size(USPEX_STRUC.POPULATION(:))
            howcome = USPEX_STRUC.POPULATION(j).howCome;
            if ~strcmp(howcome, 'keptBest') & ~strcmp(howcome, 'convexHull')
                pop_enth_full(j) = USPEX_STRUC.POPULATION(j).Enthalpies(end);
                pop_enth         = [pop_enth; USPEX_STRUC.POPULATION(j).Enthalpies(end)];
            end
            
            % Keep best can be included by uncommenting these 2 lines:
            %pop_enth_full(j) = USPEX_STRUC.POPULATION(j).Enthalpies(end);
            %pop_enth = pop_enth_full;
        end
        
        % Find enthalpies lower than target enthalpy:
        found_id_full = find(pop_enth_full <= target_enthalpy);
        found_id      = find(pop_enth      <= target_enthalpy);
        %first_id = -1;
        if ~isempty(found_id)
            first_id_full = found_id_full(1);
            first_id      = found_id(1);
            first_gen  = USPEX_STRUC.POPULATION(first_id_full).gen;
            first_enth = pop_enth(first_id);
            
            found_generations(i) = first_gen;
            found_populations(i) = first_id;
            found_enthalpies(i)  = first_enth;
            
            str = sprintf('Generation: %2i    Number: %4i    Enthalpy: %10.4f    Mat-file: %s', ...
                           first_gen,         first_id,      first_enth,         USPEX_mat{i});
            disp(str);
        end

        clear('USPEX_STRUC');
    end
    
    disp(' ');

    % Print found generations, populations, and enthalpies:
    % disp([found_populations    found_generations    found_enthalpies]);

    data_pop = [];
    data_gen = [];
    format   = [];
    for i=1:number_of_mat
        format   = [format '%4i '];
        data_pop = [data_pop found_populations(i)];
        data_gen = [data_gen found_generations(i)];
    end
    str_pop = sprintf(['Found structures numbers : ' format], data_pop);
    str_gen = sprintf(['Found generations numbers: ' format], data_gen);
    
    disp(str_pop);
    disp(str_gen)
    
    % Find not NaN values:
    positive_gen = find(found_generations > 0);
    positive_pop = find(found_populations > 0);
    
    % Calculate success rate:
    success_rate = size(positive_gen,1)/number_of_mat*100;
    
    % Calculate mean values and standard deviation:
    mean_gens = ceil(mean(found_generations(positive_gen)));
    mean_pops = ceil(mean(found_populations(positive_pop)));
    std_pops  = ceil( std(found_populations(positive_pop)));
    
    disp(' ');
    disp(['Success rate: ' num2str(success_rate) ' percent']);
    disp(['Average number of generations to get E=' num2str(target_enthalpy) ': ' num2str(mean_gens)]);
    disp(['Average number of structures  to get E=' num2str(target_enthalpy) ': ' num2str(mean_pops)]);
    disp(['Standard deviation: ' num2str(std_pops)]);
    disp(' ');
end

end
