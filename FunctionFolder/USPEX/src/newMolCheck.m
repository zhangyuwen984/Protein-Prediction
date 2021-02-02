function goodDist  = newMolCheck(MOLECULES, LATTICE, MtypeLIST,miniMatrix)

global ORG_STRUC;
% In the new version, we add a matrix to record the elmenttype for mini distance matrix
molcheck = 1;
goodDist=1;
atomType = ORG_STRUC.atomType;
STDMOL   = ORG_STRUC.STDMOL;
%status = 'checking'
%for i=1:length(minDistMOLMatrice)
%   for j=1:length(minDistMOLMatrice)
%	miniMatrix(i,j)=minDistMOLMatrice(i,j)^2;
%   end
%end
%%%%%%%%%%%check the distances between the different molecules
for i=1:length(MtypeLIST)-1
        N_atomi = length(STDMOL(MtypeLIST(i)).types);
     for j=1:N_atomi
          for k=i+1:length(MtypeLIST)
              N_atomk = length(STDMOL(MtypeLIST(k)).types);
              for l=1:N_atomk
                  dist=CalcDist(LATTICE, MOLECULES(i).MOLCOORS(j,:), MOLECULES(k).MOLCOORS(l,:));
                  type1 = STDMOL(MtypeLIST(i)).types(j);
                  type2 = STDMOL(MtypeLIST(k)).types(l);
                  if dist <= miniMatrix(type1,type2)
                  %  dist
                     goodDist=0;
                     return
                  end                    
              end       
          end
     end
end
%%%%%%%%%%%%check the distances between molecule and its periodic images
X1 = [-2:2];
Y1 = [-2:2];
Z1 = [-2:2];
[X2,Y2,Z2] = meshgrid(X1,Y1,Z1);
X3=reshape(X2,1,[]);
Y3=reshape(Y2,1,[]);
Z3=reshape(Z2,1,[]);
Matrix = [X3; Y3; Z3]';
ToDelete = all(Matrix==0, 2);
Matrix(ToDelete,:) =[]; %Remove the atoms without any contribution
%[-1 1 0; ......] except(0 0 0)

for ind=1:length(MtypeLIST)
    N_atom = length(STDMOL(MtypeLIST(ind)).types);
    for m=1:N_atom-1
        for n=m+1:N_atom
            coor1 = MOLECULES(ind).MOLCOORS(m,:);
            coor2 = MOLECULES(ind).MOLCOORS(n,:);
            tmp   = bsxfun(@plus, Matrix*LATTICE, coor2);
            dist  = min(pdist2(tmp, coor1));
            type1 = STDMOL(MtypeLIST(ind)).types(m);
            type2 = STDMOL(MtypeLIST(ind)).types(n);
            if dist <= miniMatrix(type1,type2)
               %  dist
               goodDist=0;
               return
            end
       end
    end
end
% status = 'checking done'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function dist=CalcDist(lattice, coord1, coord2)

coord1 = coord1/lattice;
coord2 = coord2/lattice;
check = coord1-coord2;
check = check - round(check);
dist = norm(check*lattice);
