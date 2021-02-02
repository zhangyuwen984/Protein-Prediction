function [X, Y] = Get_XY(numIons, comp, Energy_base, E)

comp1 = comp/numIons;
if size(numIons,2)==2
    X = comp(2)/sum(comp);
    Y = E - (comp1/sum(comp))*Energy_base;  %eV/atom
else
    if length(comp1) == 2
       X = comp1(2)/sum(comp1);
    elseif length(comp1) == 3
       tmp = comp1/sum(comp1);
       a = [-sqrt(3)/2,-1/2];  % 1 0 0
       b = [0 , 1 ];           % 0 1 0
       c = [sqrt(3)/2, -1/2];
       matrixPLOT=[a;b;c];
       X = tmp*matrixPLOT;
    end
    E = E*sum(comp)/sum(comp1);              %eV/block
    Y = E - (comp1/sum(comp1))*Energy_base;  %eV/block
end
