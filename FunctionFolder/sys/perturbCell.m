function lat = perturbCell(lattice, pertRegime)

global ORG_STRUC

% perturbs lattice to break the symmetry
% only for bulk,
pert_rate_cell = 0.5;
if size(ORG_STRUC.lattice,1) ~= 3   % Cell hasn't specified by user
    N = size(lattice, 1);
    perturbation = (rand(N,3)-0.5)*pert_rate_cell;
    lat = lattice + perturbation;
end

