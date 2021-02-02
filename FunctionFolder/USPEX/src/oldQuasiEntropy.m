function entr = oldQuasiEntropy()

% this function calculates quasi entropy that characterizes the diversity of the population
% we assume that ranking of the population is already done
% correct the case when only 1 structure is present 
% Lastly updated by Qiang Zhu (2014/02/18)

global POP_STRUC
global ORG_STRUC

N = round((length(POP_STRUC.ranking))*ORG_STRUC.bestFrac);

entr = 0;

for i = 1:N
 for j = i+1:N
   f1 = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).FINGERPRINT;
   f2 = POP_STRUC.POPULATION(POP_STRUC.ranking(j)).FINGERPRINT;
   if isempty(f1) | isempty(f2)
       cos_dist = 1;
   else
      if ORG_STRUC.varcomp == 1
        cos_dist = cosineDistance(f1, f2, 1);
      else
        cos_dist = cosineDistance(f1, f2, ORG_STRUC.weight);
      end
   end
   entr = entr + (1-cos_dist)*log(1-cos_dist);
 end
end
if N>1
   entr = 2*(-entr)/(N*N-N);
else
   entr = 0;
end
