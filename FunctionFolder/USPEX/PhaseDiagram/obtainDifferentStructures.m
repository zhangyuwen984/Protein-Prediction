function [id_list] = obtainDifferentStructures(N)
% $Rev: 789 $
% $Author: qian $
% $Date: 2014-12-03 19:45:10 +0400 (Wed, 03 Dec 2014) $

global USPEX_STRUC

id_list = [];

accepted = ones(N,1);
%--------------------------------------------------------------------------
% Fine filtering using cosine distance and rough filtering using A_order:
for i = 1 : N-1
    for j = i + 1 : N
        if SameStructure(i, j, USPEX_STRUC) == 1 || SameStructure_order(i, j, USPEX_STRUC) == 1
            accepted(j) = 0;
        end
    end
end
%--------------------------------------------------------------------------

for i = 1:N
    if accepted(i) == 1
        id_list = [id_list; i];
    end
end

end
