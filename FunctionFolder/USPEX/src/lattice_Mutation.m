function [new_Lattice,strainMatrix]= lattice_Mutation(old_lattice,mutationRate)

if nargin < 2
   %disp(['mutationRate not assgined.....']);
   mutationRate = 0.25;
end 

% lattice mutation
new_Lattice=[];
dummy = 1;
while det(new_Lattice)<0.01 | dummy 
    dummy = 0;
    strainMatrix = zeros(3);
    epsilons = randn(6,1)*mutationRate;
    strainMatrix(1,1) = 1+epsilons(1);
    strainMatrix(2,2) = 1+epsilons(2);
    strainMatrix(3,3) = 1+epsilons(3);
    strainMatrix(1,2) = epsilons(4)/2;
    strainMatrix(2,1) = epsilons(4)/2;
    strainMatrix(1,3) = epsilons(5)/2;
    strainMatrix(3,1) = epsilons(5)/2;
    strainMatrix(2,3) = epsilons(6)/2;
    strainMatrix(3,2) = epsilons(6)/2;
    new_Lattice = old_lattice*strainMatrix;
end

temp_potLat = new_Lattice;
volLat = det(temp_potLat);
% scale the lattice to the volume we asume it approximately to be
latVol = det(old_lattice);
ratio = latVol/volLat;
temp_potLat = latConverter(temp_potLat);
temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
new_Lattice = latConverter(temp_potLat);
