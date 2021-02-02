function [all_freqs, all_eigvectors, supercells_sizes] = calcSoftModes_varcomp(N_val, val1, lat, coords, numIons, maxAtoms, maxIncrease)

% calculates soft modes for different k-vectors; supercell size limited by maximum number of atoms maxAtoms 
% and maximum possible increase in any given direction maxIncrease (to avoid supercells like 1x1x100 if we start with 1 atom)
% supercells_sizes describes the supercells for corresp. freqs (aka 2x1x2, etc)

% added in 9.2.1: supercells are sorted by the length of the main diagonal, to check the smallest and bulkiest one first

global ORG_STRUC

coords0 = coords;
lat0 = lat;
numIons0 = numIons;
N = sum(numIons);
val = abs(val1);

R_val = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
 s = covalentRadius(ORG_STRUC.atomType(i));
 R_val(i) = str2num(s);
end

nMax = floor(maxAtoms/N);
nMax = min(nMax, maxIncrease);

supercells = zeros(1,3);
diagonals = 0;
% build all possible supercells
for i = 1 : nMax
 for j = 1 : nMax
  for k = 1 : nMax
    if (i*j*k*N > maxAtoms) | (i*j*k>nMax) %
      continue;
    end
    supercells = vertcat(supercells, [i j k]);
    diagonals = vertcat(diagonals, i^2 + j^2 + k^2);
  end
 end
end
supercells(1,:) = []; % remove those 0
diagonals(1) = []; % remove that 0
% sort supercells by their main diagonal
[tmp, IXd] = sort(diagonals);

all_freqs = 0;
supercells_sizes = [0 0 0];
all_eigvectors = zeros(sum(3*numIons0),1);

% calculate the frequencies/vectors for supercells
%maxCells = min(ORG_STRUC.populationSize, length(diagonals))
maxCells = length(diagonals);
for ind = 1 : maxCells
    i = supercells(IXd(ind), 1);
    j = supercells(IXd(ind), 2);
    k = supercells(IXd(ind), 3);

    %  supercell iXjXk
    kVector = zeros(1,3);
    kVector1 = zeros(1,3);
    if i > 1
      kVector(1) = 1/i;
    end
    if j > 1
      kVector(2) = 1/j;
    end
    if k > 1
      kVector(3) = 1/k;
    end
    if (i > 1) & (j > 1) & (k > 1) % 4 different sign combinations
      Nminus = 4;
      signs = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1];
    elseif (i > 1) & (j > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [1 1 0; -1 1 0];
    elseif (j > 1) & (k > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [0 1 1; 0 -1 1];
    elseif (i > 1) & (k > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [1 0 1; -1 0 1];
    else
      Nminus = 1;
      signs = [1 1 1];
    end
    for m = 1 : Nminus
      kVector1(1) = kVector(1)*signs(m,1);
      kVector1(2) = kVector(2)*signs(m,2);
      kVector1(3) = kVector(3)*signs(m,3);
      [freq, eigvector] = calcSoftModes_K(R_val, N_val, val, lat0, coords0, numIons0, kVector1);
      tmp = diag(freq);
      [sorted_freq, IX] = sort(tmp);
      for f = 1 : length(sorted_freq)
        all_freqs = vertcat(all_freqs, sorted_freq(f));
        all_eigvectors = horzcat(all_eigvectors, eigvector(:,IX(f)));
        supercells_sizes = vertcat(supercells_sizes, [i j k]);  
      end
    end
end
all_eigvectors(:,1) = []; % remove those 0
all_freqs(1,:) = []; % remove those 0
supercells_sizes(1,:) = []; % remove those 0

%experiment: let's sort
[all_freqs, IX]  = sort(all_freqs);
all_eigvectors   = all_eigvectors(:,IX);
supercells_sizes = supercells_sizes(IX,:);

%sorted_freqs = zeros(length(sorted_tmp), 4);
%for i = 1 : length(sorted_tmp)
%  sorted_freqs(i, 1) = all_freqs(IX(i), :);
%  sorted_freqs(i, 2:4) = supercells_sizes(IX(i), :);
%end

% CHANGE IT! We want to discard not the same frequencies, but same modes only!

% all_freqs = [0 0 0 0];
% tmp_f = 0;
% tmp_s = i*j*k; % supercell size
% tmp_i = 0;
% for i = 1 : size(sorted_freqs, 1)
%   if abs(sorted_freqs(i,1)) < 0.00001
%     continue;
%   end
%   cell_size = sorted_freqs(i,2)*sorted_freqs(i,3)*sorted_freqs(i,4);
%   if abs(sorted_freqs(i,1) - tmp_f) > 0.00001 % new frequency
%     if tmp_i > 0
%       all_freqs = vertcat(all_freqs, sorted_freqs(tmp_i,:)); % save the freq
%     end
%     tmp_f = sorted_freqs(i,1);
%     tmp_i = i;
%     tmp_s = cell_size;
%   elseif cell_size < tmp_s % same frequency as before, smaller cell
%     tmp_i = i;
%     tmp_s = cell_size;    
%   end
% end
% all_freqs = vertcat(all_freqs, sorted_freqs(tmp_i,:)); % save the last frequency
% all_freqs(1,:) = [];
