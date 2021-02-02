function [candidate, lattice, errorS] = Read_Stokes_output(filename, numIons)

candidate = zeros(sum(numIons),3);
lattice = [1 1 1 0 0 0]';
errorS = 0;

if exist(filename)
    try
        handle = fopen(filename);
        % fixing possible multiline minDist (only 10 numbers per output line!)
        f1 = ceil(length(numIons)^2/10);
        for LL = 1 : 6 + f1
            tmp = fgetl(handle); % repeats input file
        end
        lat_tmp = fgetl(handle); % lattice string - cell parameters:  36.84031  36.84031  36.84031  90.00000  90.00000  90.00000
        error1 = findstr(lat_tmp, 'error');
        if isempty(error1)
            lattice = sscanf(lat_tmp,'%*s %*s %g %g %g %g %g %g');
            tmp = fgetl(handle); % repeats input file
            tmp = fgetl(handle); % repeats input file
            for LL = 1 : sum(numIons)
                atom_tmp = fgetl(handle); % lattice strings in Bohr units
                atom_1 = sscanf(atom_tmp,'%*s %*s %g %g %g');
                candidate(LL,1) = atom_1(1); candidate(LL,2) = atom_1(2); candidate(LL,3) = atom_1(3);
            end
            fclose(handle);
        else
            % error1
            disp(['Error while generating crystal with symmetry ' num2str(nsym) ' and N_atoms ' num2str(numIons)]);
            fclose(handle);
            errorS = 1;
        end
    catch
        fclose(handle);
        errorS = 1;
    end
else
    disp(['ooopps, ' filename 'does not exist']);
    errorS = 1;
end
