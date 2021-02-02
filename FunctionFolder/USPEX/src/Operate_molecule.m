function R = Operate_molecule(R0, Opt, lat, P, PB)
%This routine is to apply the symmetry operations on the molcules R (in Cartesian Set)
%We need to do crystallographic operations always in fraction Set
%Also We need to consider the swap between a/b/c when it is monoclinic+orthorhombic.

%Step1: Obtain fraction axis
R0_frac = Cart2Frac(R0, lat);

%Step2: Pemutate lattice and atoms
for axis = 1:3
   R0_frac_p(:,axis) = R0_frac(:,P(axis));
   lat_p = lat(P(axis));
end

%Step3: apply crystallographic operations 
for i = 1:size(R0,1)
   R_frac_p(i,:)  = R0_frac_p(i,:)*Opt(1:3,:) + Opt(4,:);
end

%Step4: Permutate back
for axis = 1:3
    R_frac(:,axis) = R_frac_p(:,PB(axis));
end

%Step5: obtain caretesian 
R   = Frac2Cart(R_frac, lat);
