function [decomposable, samecomp, fitness] = CheckDecomposition(convex_hull, E_hull, A, Energy)
warning off
% Check the compounds for decomposition 
% -------INPUT:
% A:           input compositions
% Energy     : energy of A
% convex_hull: compositions,
% E_hull     : the energies(fitness) corrsponding to each element in the convex_hull 
% -------OUTPUT:
% decomposable : 1/0 
% samecomp     : 1/0
% fitness      : the minimum distance to the convex_hull, 

N_Comp        = size(convex_hull, 1);
N_Block       = size(convex_hull, 2);
decomposable  = 1;
samecomp      = 0;
fitness = [];
% 1: We first check if the composition is already present in the convex hull
for i = 1:N_Comp
    if sameComposition(A, convex_hull(i,:))
       fitness = Energy - E_hull(i);
       samecomp = 1;
       break;
    end
end

% 2: otherwise, we need to suggest all the possible deomposition paths, combinatorial problem
% C:  decomposition product combination of any present in the current convex_hull
% X is decomposition of A into C if CX = A has positive vector X as a solution and det(C) <> 0
% Example 1: A=[2 3]'; C=[1 0; 0 1]'; X=[2 3]'; C*X = A;
% Example 2: A=[3 2]'; C=[1 1; 0 1]'; X=[3 -1]'; C*X = A;    X must be all positive
% Example 3: A=[3 2]'; C=[1 1; 2 2]'; X would fail';         det(C) must be ~eq 0
% Example 4: A=[3 2 1]'; C=[1 0 0; 0 1 0; 1 1 0]'; failed!   det(C) must be ~eq 0

if isempty(fitness)

   Comp_List     = nchoosek([1:N_Comp],N_Block);
   N_Comp_List   = nchoosek(N_Comp,N_Block);
   C = zeros(N_Block, N_Block);
   E = zeros(N_Block, 1);
   form_Eng = [];
   N_paths = 0;
   
   for i = 1: N_Comp_List
       for j = 1: N_Block
           C(j,:) = convex_hull(Comp_List(i,j), :);
           C(j,:) = C(j,:)/sum(C(j,:));  %normalization
           E(j)   = E_hull(Comp_List(i,j));
       end
   
   if abs(det(C)) > 1.0e-8  %very important !!!
           X = A/C;
           if isempty(find(X<0)) % all X component must be positive
              N_paths = N_paths + 1;
              form_Eng(N_paths) = Energy - X*E/sum(A);  %eV/Block
           end
       end
   end

   if N_paths == 0
       disp(['Odd compositions: ' num2str(A), ' unable to find the decomposition paths']);
       fitness = 10000;
   else
       fitness = max(form_Eng);
   end
 
end

if fitness < 0
    decomposable = 0;
end
