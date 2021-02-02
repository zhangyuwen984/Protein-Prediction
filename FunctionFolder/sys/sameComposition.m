function [same] = sameComposition(c1,c2)

global ORG_STRUC

% added: Dec 26, 2010; USPEX 8.4.2
% checks whether compositions c1 and c2 are identical, 
% Here we need to put a tolerance to prevent the issue caused by a numerical error
% tolerance = 0.0001
% identity : c1 = n*c2 or c2 = n*c1 where n is integer 

L = length(c1);
tolerance = 0.0001;

ratio = 0; % compositions are identical when ratio ~= 0
for i = 1 : L
 if (abs(c1(i)) < tolerance) & (abs(c2(i)) < tolerance)
  continue;
 elseif (abs(c1(i)) < tolerance) | (abs(c2(i)) < tolerance)
  ratio = 0;
  break;
 elseif ratio == 0
  ratio = c1(i)/c2(i);
 elseif abs(c1(i)/c2(i) - ratio) > tolerance
  ratio = 0;
  break;
 end
end

if ratio == 0
 same = 0;
else
 same = 1;
end

%001mode%
% compositions are identical only if c1 = c2 (not c1 = n*c2 or c2 = n*c1)
if (ORG_STRUC.dimension==0) && (ORG_STRUC.varcomp==1) && (same == 1)
 for i = 1 : L
  if (c1(i)~=c2(i))
   same = 0;
  end
 end
end
%end of 001 mode%
