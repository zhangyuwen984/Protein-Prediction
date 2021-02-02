function [numIons, offspring, potentialLattice, fracFrac, dimension, offset, fracLattice, parents1] = heredity_final_clusterMP(parents)

global POP_STRUC
global ORG_STRUC

Np = length(parents); % number of parents
if Np > ORG_STRUC.numParents
   L = Np - ORG_STRUC.numParents;
   ind = randperm(Np);
   ind(1:L) = [];
   parents = ind;     
   Np = ORG_STRUC.numParents;
end
parentcor1 = POP_STRUC.POPULATION(parents(1)).COORDINATES(:,:);
parentcor2 = POP_STRUC.POPULATION(parents(2)).COORDINATES(:,:);
parentlat1 = POP_STRUC.POPULATION(parents(1)).LATTICE;
parentlat2 = POP_STRUC.POPULATION(parents(2)).LATTICE;
order1 = POP_STRUC.POPULATION(parents(1)).order;
order2 = POP_STRUC.POPULATION(parents(2)).order;
numIons1 = POP_STRUC.POPULATION(parents(1)).numIons;
numIons2 = POP_STRUC.POPULATION(parents(2)).numIons;
frac = 0.5;

[numIons, offspring, potentialLattice, fracFrac, dimension, offset, fracLattice] = heredity_cluster(parentcor1,parentcor2,parentlat1,parentlat2,order1,order2,numIons1,numIons2,frac);
frac = rand(1);
for i = 3 : Np
    parentcor1 = offspring;
    parentcor2 = POP_STRUC.POPULATION(parents(i)).COORDINATES(:,:);
    parentlat1 = potentialLattice;
    parentlat2 = POP_STRUC.POPULATION(parents(i)).LATTICE;
    [numIons, offspring, potentialLattice, fracFrac, dimension, offset, fracLattice] = heredity_cluster(parentcor1,parentcor2,parentlat1,parentlat2,order1,order2,frac);
end

parents1 = parents;
