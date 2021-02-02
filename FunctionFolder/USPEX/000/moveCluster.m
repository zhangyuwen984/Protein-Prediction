function [candidate] = moveCluster(lat,coordinates)
    
% this function centers the cluster in the middle of the cell
% we also put the center of mass into [0.5, 0.5, 0.5]
% first we have to center the cluster, sometimes it is divided between few neighbour cells

for j = 1 : 3
 badPositions = find(coordinates(:,j) > 0.9);
 i=0;
 while ~isempty(badPositions)
   coordinates(:,j) = coordinates(:,j) + (1-coordinates(badPositions(1),j)) + 0.01;
   coordinates = coordinates - floor(coordinates);
   badPositions = find(coordinates(:,j) > 0.9);
   i = i + 1;
   if i > 10
      %disp('we have a problem here, cluster atoms are all over the cell');
      %statusUPS = lat
      %statusUPS = coordinates
     break   % we have a huge problem here - cluster atoms are all over the cell
   end
 end
end
% center and position cluster according to the main inertia axes
AbsoluteCoord = coordinates*lat;
AbsoluteCoord = bsxfun(@minus, AbsoluteCoord, mean(AbsoluteCoord)); % Matrix-Vector;
[a,b] = PrincipleAxis(AbsoluteCoord);  %find the principle roation axis
candidate = bsxfun(@plus, AbsoluteCoord*a/lat, [0.5, 0.5, 0.5]); % move center of mass to [0,5;0.5;0.5]
