function [type, MtypeLIST, numIons]=GetPOP_MOL(numMols)

global ORG_STRUC

  MtypeLIST = [];
  for ind = 1: length(numMols)
      MtypeLIST = cat(1, MtypeLIST, ind*ones(numMols(ind),1));
  end
  numIons = [];
  type = [];
   
   for num = 1:length(MtypeLIST)
       if length(ORG_STRUC.atomType)==1
           type = cat(1,type,ORG_STRUC.atomType(ORG_STRUC.STDMOL(MtypeLIST(num)).types));
       else
           type = cat(1,type,ORG_STRUC.atomType(ORG_STRUC.STDMOL(MtypeLIST(num)).types)');
       end
   end 
   for numIND = 1: length(ORG_STRUC.atomType)
       numIons(1,numIND)=length(find(type==ORG_STRUC.atomType(numIND)));
   end 

