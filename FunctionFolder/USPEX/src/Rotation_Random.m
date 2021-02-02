function [goodSon,rotNow] = Rotation_Random(father, MtypeLIST)

global ORG_STRUC
warning off all
numMols = length(father);
rotNow  = RandInt(1,1,[1,sum(numMols)]);
STDMOL = ORG_STRUC.STDMOL;
goodSon = father;

for ind = 1: rotNow
    Coords = [];
    molRot = RandInt(1,1,[1,length(father)]);
    Coords = father(molRot).MOLCOORS;    
    num_optFlags = STDMOL(MtypeLIST(molRot)).num_optFlags;

    if(size(Coords,1)>1)
       if rand>0.8 %flip over
          Coords = bsxfun(@minus, 2*mean(Coords), Coords);  %Matrix-vector
       end

       format       = STDMOL(MtypeLIST(molRot)).format;
       [Zmatrix,MOLCOORS] = RotInertia(Coords,format, num_optFlags, MtypeLIST(molRot));
       goodSon(molRot).ZMATRIX  = Zmatrix;
       goodSon(molRot).MOLCOORS = MOLCOORS;
   end
end

