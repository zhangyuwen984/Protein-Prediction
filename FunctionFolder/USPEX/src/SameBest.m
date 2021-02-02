function same_in=SameBest()
%This program is to find the best structures has not changed for how many generations 
%Last updated by Qiang Zhu (2013/10/24)
global USPEX_STRUC
global ORG_STRUC

N_gen  = length(USPEX_STRUC.GENERATION);
N_stop = ORG_STRUC.stopCrit;
weight = ORG_STRUC.weight;
toleranceF = 0.01;

same_in = 0; 
ID1 = USPEX_STRUC.GENERATION(N_gen).BestID;

for i = N_gen - N_stop + 1 : N_gen - 1
    ID2 = USPEX_STRUC.GENERATION(i).BestID;
    if SameStructure(ID1, ID2, USPEX_STRUC)
      same_in = same_in + 1; 
    end
end
