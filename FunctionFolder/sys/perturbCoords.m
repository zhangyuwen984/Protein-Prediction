function Coords = perturbCoords(Coordinates, lattice, pertRegime)

% USPEX Version 7.2.6
% perturbs coordinates to break the symmetry
% pertRegime = 1 for atoms and clusters; pertRegime = 2 for molecules (molecules are moved as a whole); 
% pertRegime = 3 for surfaces (no z-direction movement, only non bulk part is moved)

pert_rate = 0.05; % perturbation rate in A, coords1 = coords0+-rate

N = size(Coordinates, 1);
perturbation = 2*(rand(N,3)-0.5)*pert_rate;
if pertRegime == 2   % molecules
  perturbation(:,1) = perturbation(1,1);
  perturbation(:,2) = perturbation(1,2);
  perturbation(:,3) = perturbation(1,3);
elseif pertRegime == 3    % surfaces
  perturbation(:,3) = 0;
end

perturbation = perturbation/lattice; % absolute [-0.05A;0.05A] shift goes into relative

Coords = Coordinates + perturbation;
if pertRegime ~= 2
 Coords = Coords - floor(Coords);
end