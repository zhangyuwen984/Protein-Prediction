function PreCheck_310()

global ORG_STRUC

% For fixed-cell calculations USPEX seems to have bug with non-primitive cells. % It chooses some funny cells in those cases. 
% Artem suggest to allow only P-groups for fixed-cell calculations.

%if (ORG_STRUC.nsymN(1,1) == 0) && (ORG_STRUC.constLattice)
%    for i = 1:230
%       sGroup = spaceGroups(nsym); % space group's standard symbol
%       if sGroup(1) ~= 'P'
%          ORG_STRUC.nsym(i) = 0;
%       end
%    end
%end
