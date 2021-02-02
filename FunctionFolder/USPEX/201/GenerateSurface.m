function GenerateSurface(whichInd)

global ORG_STRUC
global POP_STRUC 

CellList = findcell(ORG_STRUC.reconstruct);
Step     = POP_STRUC.POPULATION(whichInd).Step -1;
Vacuum   = ORG_STRUC.vacuumSize(Step);
  
ID = ceil(rand()*size(CellList,1));
cell=CellList(ID,:);
if ORG_STRUC.varcomp == 1
    ANS  = Random_Init_201(whichInd, Vacuum, cell);
else
    ANS  = Random_Init_200(whichInd, Vacuum, cell);
end
POP_STRUC.POPULATION(whichInd).COORDINATES         = ANS.candidate;
POP_STRUC.POPULATION(whichInd).numIons             = ANS.numIons;
POP_STRUC.POPULATION(whichInd).LATTICE             = ANS.lat;
POP_STRUC.POPULATION(whichInd).typesAList          = ANS.typesAList;
POP_STRUC.POPULATION(whichInd).chanAList           = ANS.chanAList;
POP_STRUC.POPULATION(whichInd).Surface_LATTICE     = ANS.sur_lat;
POP_STRUC.POPULATION(whichInd).Surface_COORDINATES = ANS.sur_candidate;
POP_STRUC.POPULATION(whichInd).Surface_numIons     = ANS.sur_numIons;
POP_STRUC.POPULATION(whichInd).Bulk_LATTICE        = ANS.bulk_lat;
POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES    = ANS.bulk_pos;
POP_STRUC.POPULATION(whichInd).Bulk_typesAList     = ANS.bulk_atyp;
POP_STRUC.POPULATION(whichInd).Bulk_numIons        = ANS.bulk_numIons;
POP_STRUC.POPULATION(whichInd).cell=cell;
POP_STRUC.POPULATION(whichInd).howCome = '  Random  ';
