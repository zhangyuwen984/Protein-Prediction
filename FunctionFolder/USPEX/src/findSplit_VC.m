function findSplit_VC(N, L, R1, R2, splitting)

global ORG_STRUC

% this recursive function generates all possible variations
% of atom combinations with sum between R1 and R2
% example : R1=2, R2=3 - N2, N3, O2, O3, NO2, N2O, NO
% example2 : R1=2, R2=3 - N2, N3, O2, O3, NO2, N2O, NO
% N = number of atom types, L - how many atom numbers are already determined

% USPEX 8.4.2 change : splitting determines blocks instead of number of ions

% S = sum(splitting(1:L));
S = sum(splitting(1:L)*ORG_STRUC.numIons(1:L,1:end)); % atoms = blockN*blocks

if N == L+1 % absolute recursion stop
    i1 = max(0, ceil((R1-S)/sum(ORG_STRUC.numIons(N,1:end))));
    if (S == 0) & (i1 == 0)
        i1 = 1;
    end
    i2 = floor((R2-S)/sum(ORG_STRUC.numIons(N,1:end)));
    for i = i1:i2
        splitting(N) = i;
        ORG_STRUC.splitN = ORG_STRUC.splitN + 1;
        if ORG_STRUC.splitN == 1
            ORG_STRUC.firstGeneSplit = splitting;
        else
            ORG_STRUC.firstGeneSplit = vertcat(ORG_STRUC.firstGeneSplit, splitting);
        end
    end
else
    % i1 = 0;
    i2 = floor((R2-S)/sum(ORG_STRUC.numIons(L+1,:)));
    for i = 0:i2
        splitting(L+1) = i;
        findSplit_VC(N, L+1, R1, R2, splitting) % recursion
    end
end;

ORG_STRUC.firstGeneSplitAll = ORG_STRUC.firstGeneSplit;