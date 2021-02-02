function POOL2MOL()
global POOL_STRUC
global ORG_STRUC

for pop=1:length(POOL_STRUC.POPULATION)

    numMols    = POOL_STRUC.POPULATION(pop).numMols;
    [typesAList, MtypeLIST, numIons]=GetPOP_MOL(numMols);
    lattice    = POOL_STRUC.POPULATION(pop).LATTICE;
    COORDINATES= POOL_STRUC.POPULATION(pop).COORDINATES;
    POOL_STRUC.POPULATION(pop).MtypeLIST = MtypeLIST;
    POOL_STRUC.POPULATION(pop).typesAList= typesAList;

    item = zeros(1,length(ORG_STRUC.atomType));
    newCoords = zeros(0,3);
    for mol = 1: length(MtypeLIST)
        for atom = 1: length(ORG_STRUC.STDMOL(MtypeLIST(mol)).types)
            ind = ORG_STRUC.STDMOL(MtypeLIST(mol)).types(atom);
            item(ind)=item(ind)+1;
            if ind > 1
               s = sum(numIons(1:ind-1)) + item(ind);
            else
               s = item(1);
            end
            newCoords = cat(1,newCoords, COORDINATES(s,:));
        end
    end

   for ind = 1: sum(numMols)
       amount = length(ORG_STRUC.STDMOL(MtypeLIST(ind)).types); %todo
       MOLCOORS= newCoords(1:amount,:)*lattice;
       format = ORG_STRUC.STDMOL(MtypeLIST(ind)).format;
       MOLCOORS= Movecoor(MOLCOORS, lattice, format);
       POOL_STRUC.POPULATION(pop).MOLECULES(ind).MOLCOORS= MOLCOORS;
       POOL_STRUC.POPULATION(pop).MOLECULES(ind).ZMATRIX = real(NEW_coord2Zmatrix(MOLCOORS, format));
       if length(ORG_STRUC.STDMOL(MtypeLIST(ind)).types) ==1
          POOL_STRUC.POPULATION(pop).MOLECULES(ind).MOLCENTER = MOLCOORS;
       else
          POOL_STRUC.POPULATION(pop).MOLECULES(ind).MOLCENTER = mean(MOLCOORS);
       end
       newCoords(1:amount,:)=[];
       MOLCOORS = [];
   end
end


function newCoor = Movecoor(coords,lat,format)

newCoor = coords;
if(size(coords,1)>1)
 for ind = 2:size(coords,1)
   check = newCoor(ind,:)-newCoor(format(ind,1),:);
   mindist = 10;
   for x = 1:3
      for y=1:3
         for z=1:3
             dist = sqrt(sum((check+[x-2,y-2,z-2]*lat).^2));
             if dist < mindist
                mindist = dist;
                newCoor(ind,:) = coords(ind,:) + [x-2,y-2,z-2]*lat;
             end
         end
      end
   end
 end
end

