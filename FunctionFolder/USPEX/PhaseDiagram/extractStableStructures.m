function extractStableStructures(ids)
% $Rev: 1120 $
% $Author: mrakitin $
% $Date: 2015-08-20 21:51:42 +0400 (Thu, 20 Aug 2015) $

global ORG_STRUC
global USPEX_STRUC

% Save found structures to stablePOSCARS_pressure:
[nothing, nothing] = unix('rm -f POSCAR stablePOSCARS_pressure');

for i=1:length(ids)
    Write_POSCAR(ORG_STRUC.atomType, ids(i),             ...
                 USPEX_STRUC.POPULATION(ids(i)).symg,    ...
                 USPEX_STRUC.POPULATION(ids(i)).numIons, ...
                 USPEX_STRUC.POPULATION(ids(i)).LATTICE, ...
                 USPEX_STRUC.POPULATION(ids(i)).COORDINATES);

    [nothing, nothing] = unix('cat POSCAR >> stablePOSCARS_pressure');
end

[nothing, nothing] = unix('rm -f POSCAR');

end
