function [lat1, coor1]=reduce2D(lat, coor, thickness)
lat1 = lat;
lat1(3,3) = thickness;
lat1(1,3) = 0;
lat1(2,3) = 0;
lat1(3,1) = 0;
lat1(3,2) = 0;
coor1 = coor*lat/lat1;
aveh = sum(coor1(:,3))/size(coor1,1);
shift = aveh - 0.5;
for i = 1:size(coor1,1)
   coor1(i,3) = coor1(i,3) - shift;
end
coor1 = coor1 - floor(coor1);

