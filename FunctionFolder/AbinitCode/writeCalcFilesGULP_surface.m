function writeCalcFilesGULP_surface(Ind_No)

% USPEX Version 8.1.1
% Change: added surface support

global POP_STRUC
global ORG_STRUC

LATTICE         = POP_STRUC.POPULATION(Ind_No).LATTICE;
COORDINATES     = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons         = POP_STRUC.POPULATION(Ind_No).numIons;
Surface_numIons = POP_STRUC.POPULATION(Ind_No).Surface_numIons;
Bulk_numIons    = POP_STRUC.POPULATION(Ind_No).Bulk_numIons;

if sum(numIons)~=sum(Surface_numIons+Bulk_numIons)
  disp('EEEEEEEEEEEEE writeCalc')
  quit
end

COORDINATES = COORDINATES*LATTICE;

[nothing, nothing] = unix(['echo svectors >> input']);

[nothing, nothing] = unix(['echo ' num2str(Lattice(1,1)) '  ' num2str(LATTICE(1,2)) ' >> input']);
[nothing, nothing] = unix(['echo ' num2str(Lattice(2,1)) '  ' num2str(LATTICE(2,2)) ' >> input']);

[nothing, nothing] = unix(['echo ''cartesian region 1'' >> input']);
item=0;

for i=1:size(COORDINATES,1)
  if POP_STRUC.POPULATION(Ind_No).chanAList(i)==1
     [nothing, nothing] = unix(['echo ' num2str(POP_STRUC.POPULATION(Ind_No).typesAList(i)) ' ' num2str(COORDINATES(i,:)) '  >> input']);
     item=item+1;
     POP_STRUC.POPULATION(Ind_No).Index(item)=i;
  end
end

[nothing, nothing] = unix(['echo ''cartesian region 2'' >> input']);
for i=1:size(COORDINATES,1)
   if POP_STRUC.POPULATION(Ind_No).chanAList(i)<1
      [nothing, nothing] = unix(['echo ' num2str(POP_STRUC.POPULATION(Ind_No).typesAList(i)) ' ' num2str(COORDINATES(i,:)) '  >> input']);
      item=item+1;
      POP_STRUC.POPULATION(Ind_No).Index(item)=i;
   end
end
