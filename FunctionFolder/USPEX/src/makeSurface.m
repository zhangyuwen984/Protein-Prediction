function [lat, coordinates, numIons, typesAList, chanAlist] = makeSurface(sur_lattice, sur_coordinate, sur_numIons, bulk_lattice, bulk_pos, bulk_atyp, bulk_numIons, vacuumSize)
global ORG_STRUC 
  bulk_lat = latConverter(bulk_lattice);
  sur_lat  = latConverter(sur_lattice);

  tempsurf(1:size(sur_coordinate,1)) = 1;
  tempbulk(1:size(bulk_pos,1)) = 0;

  lat = sur_lat;
  lat(3) = sur_lat(3)+bulk_lat(3)+vacuumSize; %total lat

  if size(sur_coordinate,1) ~=0
     sur_coordinate(:,3) = (bulk_lat(3)+sur_coordinate(:,3)*sur_lat(3))/lat(3); 
  end

  bulk_pos(:,3) = bulk_pos(:,3)*bulk_lat(3)/lat(3);


  [lattice] = latConverter(lat);
  lat = lattice;                       %%total lattice
  coord = cat(1,sur_coordinate,bulk_pos);  %%total coordinates
  alist = cat(2,tempsurf, tempbulk);   %%total list
  numIons = sur_numIons+bulk_numIons; %%total Ions

sur_list=[];
for i = 1:size(sur_numIons,2)
  item = length(sur_list);
  if sur_numIons(i)>0
     sur_list(item+1 : item + sur_numIons(i)) = ORG_STRUC.atomType(i);
  end
end

typesAList = cat(2,sur_list,bulk_atyp');
atomlists = typesAList;
k = 0;
for j = 1:length(ORG_STRUC.atomType)
 for i = 1:size(coord,1)
  if typesAList(i) == ORG_STRUC.atomType(j)
   k = k+1;
   coordinates(k,:) = coord(i,:);
   chanAlist(k) = alist(i);
   atomlists(k) = typesAList(i);
  end
 end
end
typesAList = atomlists;
