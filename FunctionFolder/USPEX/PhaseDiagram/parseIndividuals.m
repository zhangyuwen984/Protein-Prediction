function [comp, enth, volume] = parseIndividuals(id_list)
% $Rev: 685 $
% $Author: maxim $
% $Date: 2014-10-24 10:30:20 +0400 (Fri, 24 Oct 2014) $

global ORG_STRUC
global USPEX_STRUC

comp   = [];
enth   = [];
volume = [];
for i=1:length(id_list)
    if ORG_STRUC.molecule == 1
        comp   = [comp; USPEX_STRUC.POPULATION(i).numMols];
    else
        comp   = [comp; USPEX_STRUC.POPULATION(i).numIons];
    end
    enth   = [enth  ; USPEX_STRUC.POPULATION(i).Enthalpies(end)];
    volume = [volume; USPEX_STRUC.POPULATION(i).Vol];
end

end
