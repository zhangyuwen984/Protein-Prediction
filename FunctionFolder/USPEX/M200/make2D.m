function [lat, coordinates] = make2D(lattice1, coordinate1, vacuumSize)
  lat  = lattice1;
  lat(3,3) = lat(3,3)+vacuumSize; %total lat
  lat(1,3) = 0;
  lat(2,3) = 0;
  lat(3,1) = 0;
  lat(3,2) = 0;
  coordinates = coordinate1;
  coordinates = coordinates - floor(coordinates);
  for i = 1:size(coordinate1,1)
     coordinates(i,3) = coordinates(i,3) - 0.5;  %to move to [-0.5:0.5]
     coordinates(i,3) = 0.5 + coordinates(i,3)*lattice1(3,3)/lat(3,3);
  end

