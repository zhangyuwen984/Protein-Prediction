function MakeupSurface(whichInd,lattice,coor,type,surnumIons)

%%%%%%%%%This routine is used to compensate the atom if some add atoms are removed from the network during relaxation%%%%%%%%%%%
global ORG_STRUC
global POP_STRUC 

Step = POP_STRUC.POPULATION(whichInd).Step;
numatom = ORG_STRUC.numIons*prod(POP_STRUC.POPULATION(whichInd).cell) - surnumIons; 
lattice = POP_STRUC.POPULATION(whichInd).Surface_LATTICE;
bulk_lattice=POP_STRUC.POPULATION(whichInd).Bulk_LATTICE;
bulk_pos=POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES;
bulk_atyp=POP_STRUC.POPULATION(whichInd).Bulk_typesAList;
bulk_numIons=POP_STRUC.POPULATION(whichInd).Bulk_numIons;          
[lat,candidate,numIons,typesAList,chanAList] = makeSurface(lattice,coor,surnumIons,bulk_lattice,bulk_pos,bulk_atyp,bulk_numIons,ORG_STRUC.vacuumSize(Step)); 
Occupied=zeros(1,length(ORG_STRUC.numIons));
adatom_type = [];
adatom_coord = [];

 for i=1:sum(numatom)
     goodpick = 0;
     while ~goodpick
%%%%%try to pick the atom to be added randomly
       t=RandInt(1,1,[1,length(ORG_STRUC.atomType)]);
       if Occupied(t) < numatom(t)+0.5
          goodpick=1;
          [x,y]=choose_docking([adatom_coord;candidate],lat,[adatom_type;typesAList']',ORG_STRUC.atomType(t));
          zmin = bulk_lattice(3,3)/lat(3,3);
          z = surface_adatom_docking(lat,zmin,[adatom_coord;candidate],[adatom_type;typesAList'],[x,y],ORG_STRUC.atomType(t));
          adatom_coord=[adatom_coord;[x,y,z]];
          adatom_type = [adatom_type;ORG_STRUC.atomType(t)];
          Occupied(t)=Occupied(t)+1;
       end
     end
 end
 adatom_type = [adatom_type;type'];
 adatom_coord = [adatom_coord;coor];
 for i=1:sum(numatom)
     adatom_coord(i,3) = (adatom_coord(i,3)*lat(3,3)-bulk_lattice(3,3))/lattice(3,3);
 end     
     newCoords = zeros(0,3);
 for ind = 1: length(ORG_STRUC.atomType)
     s = find(adatom_type == ORG_STRUC.atomType(ind));
     newCoords = cat(1,newCoords, adatom_coord(s,:));
 end

 surnumIons1 = ORG_STRUC.numIons*prod(POP_STRUC.POPULATION(whichInd).cell);
 [lat,candidate,numIons,typesAList,chanAList] = makeSurface(lattice,newCoords,surnumIons1,bulk_lattice,bulk_pos,bulk_atyp,bulk_numIons,ORG_STRUC.vacuumSize(Step));

 POP_STRUC.POPULATION(whichInd).Surface_COORDINATES = newCoords;
 POP_STRUC.POPULATION(whichInd).Surface_LATTICE = lattice;
 POP_STRUC.POPULATION(whichInd).Surface_numIons = surnumIons1;
 POP_STRUC.POPULATION(whichInd).COORDINATES = candidate;
 POP_STRUC.POPULATION(whichInd).numIons = numIons;
 POP_STRUC.POPULATION(whichInd).LATTICE = lat;
 POP_STRUC.POPULATION(whichInd).typesAList = typesAList;
 POP_STRUC.POPULATION(whichInd).chanAList=chanAList;
