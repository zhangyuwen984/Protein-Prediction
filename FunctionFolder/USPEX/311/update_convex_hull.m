function [fitness] = update_convex_hull()

% if N_T == 1 - it will fail!

global POP_STRUC
global ORG_STRUC

N_T = size(ORG_STRUC.numIons,1);
N_type = size(ORG_STRUC.numIons,2);
N_P = length(POP_STRUC.POPULATION);
numIons = zeros(N_P, N_T);
for it=1:N_P
   for it2=1:N_type
         numIons(it,it2) = POP_STRUC.POPULATION(it).numMols(it2);
   end
end
fitness = zeros(1, N_P);

% it turns out that gulp can receive worse values for convex hill points after recalculation
% therefore we can't use old convex_hull and have to build new one from scratch it seems

%if POP_STRUC.generation == 1
 POP_STRUC.convex_hull = zeros(N_T, N_T + 2);
 for type_loop = 1 : N_T
   POP_STRUC.convex_hull(type_loop, N_T + 1) = 1000000;
 end
%end
ch_add = zeros(1, N_T + 2);
finished = zeros(1, N_P);


% find the elementary substances
for type_loop = 1 : N_T
  for it = 1 : N_P
    if POP_STRUC.POPULATION(it).numBlocks(type_loop) == sum(POP_STRUC.POPULATION(it).numBlocks)
       fitn = POP_STRUC.POPULATION(it).Enthalpies(end)/sum(numIons(it, :));
       finished(it) = 1;
       if fitn <= POP_STRUC.convex_hull(type_loop, N_T + 1) + 0.0000001
          POP_STRUC.convex_hull(type_loop, 1:N_T) = POP_STRUC.POPULATION(it).numBlocks;
          POP_STRUC.convex_hull(type_loop, N_T + 1) = fitn;
          POP_STRUC.convex_hull(type_loop, N_T + 2) = it;
       end
    end
  end
end

tic

% check already known compositions
for c_loop = 1 : size(POP_STRUC.convex_hull,1)
  for it = 1 : N_P
    same = sameComposition(POP_STRUC.convex_hull(c_loop, 1:N_T), POP_STRUC.POPULATION(it).numBlocks);
    if same == 0
     continue;
    else 
     fitn = POP_STRUC.POPULATION(it).Enthalpies(end)/sum(numIons(it, :));
     finished(it) = 1;
     if fitn - POP_STRUC.convex_hull(c_loop, N_T + 1) < 0.000001
       POP_STRUC.convex_hull(c_loop, 1:N_T) = POP_STRUC.POPULATION(it).numBlocks;
       POP_STRUC.convex_hull(c_loop, N_T+1) = fitn;
       POP_STRUC.convex_hull(c_loop, N_T+2) = it;
     end
    end
  end
end

C = zeros(N_T, N_T);

% Here we check the compounds for decomposition, this is done via binary string where 1 means that structure is in the basis for decomposition 
% those elements build matrix C with their atom coefficients; X is decomposition of A into C if CX = A has positive vector X as a solution and det(C) <> 0
% First check decomposition for single elements

for it = 1 : N_P

 if finished(it)
  continue;
 end
 A = POP_STRUC.POPULATION(it).numBlocks;
 A = A'; % make it a column
 compos = zeros(1, size(POP_STRUC.convex_hull,1));
 for i = 1 : N_T
  compos(i) = 1;
  C(i,:) = POP_STRUC.convex_hull(i, 1:N_T);
 end
 C = C';
 fitn = POP_STRUC.POPULATION(it).Enthalpies(end)/sum(numIons(it, :));
 decomposable = 0;

 while 1 % convex_hull cycle
  
  if det(C) ~= 0
    %X = C\A;
    X = pinv(C)*A;
     if sum(X<-0.001)==0 % all X component must be positive
      f = 0;
      k = 1;
      for i = 1 : size(POP_STRUC.convex_hull,1) 
       if compos(i) 
         f = f + X(k)*POP_STRUC.convex_hull(i, N_T + 1)*sum((POP_STRUC.convex_hull(i, 1:N_T)*ORG_STRUC.numIons));
         k = k + 1;
       end
      end
      f = f/sum(numIons(it, :));
      if fitn > f % fitnesswise can be decomposed
        decomposable = 1;
        break;
      end
    end
  end

  flag = 1;
  for i = size(POP_STRUC.convex_hull,1) - N_T + 1 : size(POP_STRUC.convex_hull,1)
   if compos(i) == 0
     flag = 0;
   end
  end
  if flag
   break;
  end
% update compos
  k = size(POP_STRUC.convex_hull,1) - 1;
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
  compos(k+s+1 : size(POP_STRUC.convex_hull,1)) = 0;
  k = 1;
  for i = 1 : size(POP_STRUC.convex_hull,1)
   if compos(i)
    C(k,:) = POP_STRUC.convex_hull(i, 1:N_T);
    k = k + 1;
   end
  end
  C = C';
 end

 if decomposable == 0
   ch_add(1:N_T) = POP_STRUC.POPULATION(it).numBlocks;
   ch_add(N_T+1) = fitn;
   ch_add(N_T+2) = it;
   POP_STRUC.convex_hull = vertcat(POP_STRUC.convex_hull, ch_add);
  % check identical compositions
  for c_loop = 1 : size(POP_STRUC.convex_hull,1)
   for it1 = 1 : N_P
    same = sameComposition(POP_STRUC.convex_hull(c_loop, 1:N_T), POP_STRUC.POPULATION(it1).numBlocks);
    if same == 0
     continue;
    else 
     fitn = POP_STRUC.POPULATION(it1).Enthalpies(end)/sum(numIons(it1, :));
     finished(it1) = 1;
     if fitn - POP_STRUC.convex_hull(c_loop, N_T + 1) < 0.000001
       POP_STRUC.convex_hull(c_loop, 1:N_T) = POP_STRUC.POPULATION(it1).numBlocks;
       POP_STRUC.convex_hull(c_loop, N_T+1) = fitn;
       POP_STRUC.convex_hull(c_loop, N_T+2) = it1;
     end
    end
   end
  end
 end
 finished(it) = 1;

end

t = toc;
disp(['New points added to convex hull in ' num2str(t) ' sec'])
disp(' ')

% Calculating fitness = min(fitn - any possible decomposition energy)
% To check the compounds for decomposition, we use binary string where 1 means that structure is in the basis for decomposition 
% those elements build matrix C with their atom coefficients; X is decomposition of A into C if CX = A has positive vector X as a solution and det(C) <> 0

for it = 1 : N_P

 A = POP_STRUC.POPULATION(it).numBlocks;
 A = A'; % make it a column
 compos = zeros(1, size(POP_STRUC.convex_hull,1));
 for i = 1 : N_T
  compos(i) = 1;
  C(i,:) = POP_STRUC.convex_hull(i, 1:N_T);
 end
 fitn = POP_STRUC.POPULATION(it).Enthalpies(end)/sum(numIons(it, :));
 fitness(it) = 1000000;
 split = compos;
 cc = compos;
 C = C';

 while 1 % convex_hull cycle
  
  if det(C) ~= 0
    %X = C\A;
    X = pinv(C)*A;
    if isempty(find(sign(X) == -1)) % composition can be decomposed
      f = 0;
      k = 1;
      for i = 1 : size(POP_STRUC.convex_hull,1) 
       if compos(i) 
         f = f + X(k)*POP_STRUC.convex_hull(i, N_T + 1)*sum((POP_STRUC.convex_hull(i, 1:N_T)*ORG_STRUC.numIons));
         k = k + 1;
       end
      end
      f = f/sum(numIons(it, :));
      if (f - fitn) < fitness(it)
        fitness(it) = f - fitn;
        split = X;
        cc = compos;
      end
    end
  end

  flag = 1;
  for i = size(POP_STRUC.convex_hull,1) - N_T + 1 : size(POP_STRUC.convex_hull,1)
   if compos(i) == 0
     flag = 0;
   end
  end
  if flag
   break;
  end
% update compos
  k = size(POP_STRUC.convex_hull,1) - 1;
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
  compos(k+s+1 : size(POP_STRUC.convex_hull,1)) = 0;
  k = 1;
  for i = 1 : size(POP_STRUC.convex_hull,1)
   if compos(i)
    C(k,:) = POP_STRUC.convex_hull(i, 1:N_T);
    k = k + 1;
   end
  end
  C = C';
 end

% remark it if you don't want 'parabolic' correction:
% k1 = 0;
% k = 1;
% for i = 1 : size(POP_STRUC.convex_hull,1) 
%  if cc(i) 
%   if split(k) > 0.00001
%     k1 = k1 + 1;
%     fitness(it) = fitness(it)/split(k);
%   end
%   k = k + 1;
%  end
% end
% fitness(it) = fitness(it)/(k1^k1);

end

% update convex hull - throw all those points out, that have fitness not equal to zero

c_loop = 1;

while 1
 if c_loop > size(POP_STRUC.convex_hull,1)
  break;
 end
 for it = 1 : N_P
    diff = sum(abs(POP_STRUC.convex_hull(c_loop, 1:N_T) - POP_STRUC.POPULATION(it).numBlocks));
    fitn = POP_STRUC.POPULATION(it).Enthalpies(end)/sum(numIons(it, :));
    if (diff == 0) & (abs(fitn - POP_STRUC.convex_hull(c_loop, N_T + 1)) < 0.0001)
     if abs(fitness(it)) > 0.0001
      POP_STRUC.convex_hull([c_loop],:) = []; 
      c_loop = c_loop - 1;
     end 
     break;
    end
 end
 c_loop = c_loop + 1;
end

t = toc;
disp(['Convex hull updated and fitness calculated in ' num2str(t) ' sec'])
disp(' ')

fitness = -1*fitness;
