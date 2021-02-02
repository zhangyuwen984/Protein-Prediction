function content = generate_poscar(lattice, coordinates, atoms, num_atoms, file, nsym)
% GENERATE_POSCAR: generate POSCAR file.

if nargin < 5
    file = 'POSCAR.vasp';
    nsym = '0';  % default value for symmetry group if not specified
end

%% Create cells for lattice and coordinates:
lat = {};
for i=1:3
    lat(end+1) = {num2str(lattice(i,:), '%12.5f')};
end

crds = {};
for i = 1:sum(num_atoms)
    crds(end+1) = {num2str(coordinates(i,:), '%12.5f')};
end


%% Create content cell variable:
content = {
    ['POSCAR structure with SG #: ' num2str(nsym)]
    '1.0'
    };

for i=1:size(lat, 2)
    content(end+1) = lat(i);
end

content(end+1) = atoms;
content(end+1) = {num2str(num_atoms)};
content(end+1) = {'Direct'};

for i=1:size(crds, 2)
    content(end+1) = crds(i);
end


%% Print data to a file:
f = fopen(file, 'w');
for i=1:size(content, 1)
    fprintf(f, '%s\n', char(content(i)));
end
fclose(f);
