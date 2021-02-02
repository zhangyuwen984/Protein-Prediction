function [lat, candidate] = makeCluster(lattice, coordinates, vacuumSize)
% $Rev: 584 $
% $Author: maxim $
% $Date: 2014-08-30 10:27:47 +0400 (Sat, 30 Aug 2014) $

% this function adds 'free space' to the cluster, evenly from all sides

lat = lattice;
coordinates = coordinates*lat; 

lat(1,1) = lat(1,1) + vacuumSize;
lat(2,2) = lat(2,2) + vacuumSize;
lat(3,3) = lat(3,3) + vacuumSize;

coordinates = coordinates + vacuumSize/2;

coordinates = coordinates/lat;

candidate = coordinates;
