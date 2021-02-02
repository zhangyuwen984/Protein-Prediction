function [newLat, newCoord] = rotateCluster(coord, lat, alpha, beta, gamma)
% USPEX 9.2.3
% rotates the cluster around axis X for alpha, around Y for beta and around Z for gamma
% then 'centers' the cluster in 0.5 0.5 0.5 and builds a cell around it
AbsCoord = coord*lat;
mass_center = zeros(1,3);

AbsCoord = bsxfun(@minus, AbsCoord, mean(AbsCoord)); % Matrix - Vector

Cnz = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1]; % rotate by Z axis 
Cnx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]; % rotate by X axis
Cny = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];     % rotate by Y axis

% change the lattice !!!!!
AbsCoord = AbsCoord*Cnz*Cnx*Cny;

newCoord = bsxfun(@plus, AbsCoord/lat, [0.5, 0.5, 0.5]); % move center of mass to [0,5;0.5;0.5]
[newLat, newCoord] = reduce_Cluster(lat, newCoord);

