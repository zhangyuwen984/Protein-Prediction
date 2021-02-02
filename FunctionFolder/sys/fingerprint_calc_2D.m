function [order, fing, atom_fing]=fingerprint_calc_2D(Ni, V, dist_matrix, typ_i, typ_j, numIons, ho, ht)

% Ni - number of atoms in the unit cell of sort that is determined by atom i
% Ncell - number of atoms in the unit cell
% Nfull - full number of atoms that we use in our analysis
% typ_i,j - types of atoms i and j (1,2,3,...,species)
% Reference: AR Oganov, M Valle, How to quantify energy landscapes, J. Chem. Phys, 104504, 2009
global ORG_STRUC

Rmax  = ORG_STRUC.RmaxFing;
sigm  = ORG_STRUC.sigmaFing;
delta = ORG_STRUC.deltaFing;
species = size(numIons,2);
Ncell = sum(numIons);
normaliser = 1;
Nfull = size(dist_matrix, 1);

if sum(Ncell)==0
      order=[];
      fing =[];
      atom_fing=[];
else
       numBins = round(Rmax/delta);
       fing = zeros(species*species, numBins);
       fing(:,1) = -1;

       order = zeros(1, Ncell);
       atom_fing = zeros(Ncell, species, numBins);
       sigm = sigm/sqrt(2*log(2));
       
       interval=zeros(2,1);
end       
       for bins = 2:numBins
        for i = 1:Ncell
            for j = 1:Nfull
               if (dist_matrix(j,i)>0) & (abs(dist_matrix(j,i)-delta*(bins-0.5))<4*sigm)
                   R0 = dist_matrix(j,i);
       
                   interval(2)=5*sign(delta*bins-R0);
                   if abs((delta*bins-R0)/sqrt(2)) <= 5*sigm
                       interval(2) = (delta*bins-R0)/(sqrt(2)*sigm);
                   end
                   interval(1)=5*sign(delta*(bins-1)-R0);
                   if abs((delta*(bins-1)-R0)/sqrt(2)) <= 5*sigm
                       interval(1) = (delta*(bins-1)-R0)/(sqrt(2)*sigm);
                   end
                   down = ho(1,i);
                   up = ht (1,i);                     
                   my_erf=approx_erf(interval, ORG_STRUC.erf_table);
                   delt=0.5*(my_erf(2)-my_erf(1));
                   atom_fing(i, typ_j(j), bins) = atom_fing(i, typ_j(j), bins) + delt/(Ni(j)*R0^2); % delt - gaussian with half-width sigma and norm 1
                  if abs(R0) <= down && abs(R0) <= up
                  fing((typ_i(i)-1)*species+typ_j(j), bins) = fing((typ_i(i)-1)*species+typ_j(j), bins) + delt/(2*R0^2); % delt - gaussian with half-width sigma and norm 1
                  
                  elseif abs(R0) < up && abs(R0) > down
                  fing((typ_i(i)-1)*species+typ_j(j), bins) = fing((typ_i(i)-1)*species+typ_j(j), bins) + delt/(R0^2+R0*down); % delt - gaussian with half-width sigma and norm 1
                  
                  elseif abs(R0) < down && abs(R0) > up
                  fing((typ_i(i)-1)*species+typ_j(j), bins) = fing((typ_i(i)-1)*species+typ_j(j), bins) + delt/(R0^2+R0*up); % delt - gaussian with half-width sigma and norm 1
                                 
                  elseif abs(R0) > down && abs(R0) > up
                  fing((typ_i(i)-1)*species+typ_j(j), bins) = fing((typ_i(i)-1)*species+typ_j(j), bins) + delt/(R0*(up+down)); % delt - gaussian with half-width sigma and norm 1
                  end
               end
            end
            atom_fing(i, :, bins) = V*atom_fing(i, :, bins)/(4*pi*delta) - normaliser;
            for j = 1:species
                weight = numIons(j)/sum(numIons); %weight factor applied by Qiang Zhu (2013/11/04)
                order(i) = order(i) + weight * delta*atom_fing(i, j, bins)*atom_fing(i, j, bins)/power(V/Ncell,1/3);  % eq(5) in CPC-2010 paper
            end
        end
       
       for i = 1:species
         for j = 1:species
             if numIons(i)>0 & numIons(j)>0
           fing((i-1)*species+j, bins) = V*fing((i-1)*species+j, bins)/(2*pi*numIons(i)*numIons(j)*delta) - normaliser;
         end;
       end;
       
       end
       order = sqrt(order);
end
