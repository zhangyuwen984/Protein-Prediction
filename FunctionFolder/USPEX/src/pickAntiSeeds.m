function pickAntiSeeds()

% USPEX Version 7.4.3
% Change: added

global ANTISEEDS
global ORG_STRUC

ANTISEEDS = struct('FINGERPRINT',{},'Sigma',{},'Max',{});
ANTISEEDS(1).FINGERPRINT = [];
ANTISEEDS(1).Sigma = [];
ANTISEEDS(1).Max = [];

cd AntiSeeds

disp('Read AntiSeeds ...')

try
    % reading POSCARS

    [fid,message] = fopen('POSCARS');
    loops_seed = 0;
    while 1

            loops_seed = loops_seed + 1;
        %%%=================================%%%
        %%% Read the AntiSeeds POSCAR file  %%%
        %%%=================================%%%

            tmp = fgetl(fid); % system description

            scale_factor = fgetl(fid);           % 1 = numbers in angstrems
            optlattice = fscanf(fid,'%g',[3,3]); % lattice vectors
            optlattice = optlattice'*str2num(scale_factor);
            tmp = fgetl(fid);                    % we don't need this line

            atomType = fgetl(fid);   % atomic element symbol
            ntyp = fgetl(fid);       % types of atoms aka atomic numbers
            ntyp = str2num(ntyp);

            natom = sum(ntyp);       % number of atoms
            tmp = fgetl(fid);        % Direct/Cart line
            candidate_pop = fscanf(fid,'%g',[3,natom]);
            candidate_pop = candidate_pop';
            tmp = fgetl(fid);


        %%%==========================================%%%
        %%% Calculate the FingerPrint with AntiSeeds %%%
        %%%==========================================%%%

            numIons = ntyp;

            [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(optlattice, candidate_pop, numIons, 1);
            [order, fing, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);

            ANTISEEDS(loops_seed).FINGERPRINT = fing;

            ANTISEEDS(loops_seed).Sigma = ORG_STRUC.antiSeedsSigma;
            ANTISEEDS(loops_seed).Max   = ORG_STRUC.antiSeedsMax;

            disp([' --> AntiSeed ' num2str(loops_seed) ' added ...'])
    end
    status = fclose(fid);

catch
    %%% Do nothing here.
end

fclose('all');
cd (ORG_STRUC.homePath)

% save the AntiSeeds
safesave ('ANTISEEDS.mat', ANTISEEDS);