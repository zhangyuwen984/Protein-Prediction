function [fitness, convex_hull0] = phase_update_hull(N_Block, N_Point, numBlocks, Energy, numIons)

% if N_T == 1 - it will fail!
% Completely rewrite this code
% 

% All input variables are passed to the function arguments.

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- Convex hull calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it turns out that gulp can receive worse values for convex hill points after recalculation
% therefore we can't use old convex_hull and have to build new one every generation

fitness = zeros(1, N_Point);
convex_hull = 1000000*ones(N_Block, N_Block + 2);
ch_add      = zeros(1, N_Block + 2);
finished    = zeros(1, N_Point);

% find the elementary substances (to define the boudary of convex hull)
% The basis set always starts from [1 0 0 ...] [0 1 0 .....], ...
% The dimension depends on N_Block
% Other compositions thus can be reprenseted by those basis set
% One can thus use the following to calulate other convex hull problem

for i = 1 : N_Point
    for j = 1: N_Block
    if numBlocks(i,j) == sum(numBlocks(i,:))
       finished(i) = 1;
       if Energy(i) <= convex_hull(j, N_Block+1) + 0.0000001
          convex_hull(j, 1 : N_Block)   = numBlocks(i,:);  %compositions, 
          convex_hull(j, N_Block + 1)   = Energy(i);       %Energies
          convex_hull(j, N_Block + 2)   = i;               %Index, just for convenience
       end
    end
  end
end

% Calculate other compositions

for i = 1 : N_Point
    if finished(i)
       continue;
    end

    [decomposable, samecomp, fitness(i)] = CheckDecomposition(convex_hull(:, 1:N_Block), ...
                             convex_hull(:, N_Block+1), numBlocks(i,:), Energy(i));

    if decomposable == 0 %If the new structure is stable, we need to do the following
       ch_add  = [numBlocks(i,:) Energy(i) i];
       if samecomp == 1  %1, if the composition already exists in the convex hull, REPLACE!
          for j = 1:size(convex_hull, 1)
              if sameComposition(numBlocks(i,:), convex_hull(j,:))
                 convex_hull(j,:) = ch_add;
                 break
              end
          end
       else              %2, if new composition, ADD!
           convex_hull   = [convex_hull; ch_add];
       end

       % check if other compositions in the convex hull will be impacted 
       % Here we check each composition one by one
       % If valid, we include it
       convex_hull1 = convex_hull;
       convex_hull(N_Block+1:end, :) = [];
       for j = N_Block+1 : size(convex_hull1, 1)
           tmp_hull = convex_hull1;
           tmp_hull(j,:) = [];
           [tmp1, tmp2] = CheckDecomposition(tmp_hull(:, 1:N_Block), ...
           tmp_hull(:, N_Block+1), convex_hull1(j, 1:N_Block), convex_hull1(j, 1+N_Block));
           if tmp1 == 0  % keep the survived composition
              convex_hull = [convex_hull; convex_hull1(j,:)];
           end
       end
    end
end

% Now we recalculate the fitness
for i = 1 : N_Point
    [decomposable, samecomp, fitness(i)] = CheckDecomposition(convex_hull(:, 1:N_Block), ...
                             convex_hull(:, N_Block+1), numBlocks(i,:), Energy(i));
end

% In USPEX, we output the convex hull by elemental composition, but not block composition
convex_hull0 = Transform_Convex_Hull(convex_hull, numIons, 0);

t = toc;
%disp(['Calculating fitness in ' num2str(t) ' sec'])
%disp(' ')
