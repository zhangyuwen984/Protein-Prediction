function [T_Lat, T_coor, T_atyp, chanAList] = makeGB(numIons, GB_lat, GB_coor, GB_atyp, bulk_lat, bulk_pos, bulk_atyp, vacuum)

global ORG_STRUC 

GB_num = length(GB_atyp);
bulk_num = length(bulk_atyp');
T_coor = [];
chanAList = [];
lat1 = latConverter(GB_lat);
lat2 = latConverter(bulk_lat);
lat = lat1;
lat(3) = vacuum+lat2(3); %total lat
T_Lat = latConverter(lat);

minH = (lat2(3)-lat1(3))/2;

bulk_pos(:,3) = bulk_pos(:,3)*lat2(3)/lat(3);

for i = 1:size(GB_coor)
   GB_coor(i,3)  = GB_coor(i,3)*lat1(3) + minH;
   GB_coor(i,3)  = GB_coor(i,3)/lat(3);
end

coord = [GB_coor; bulk_pos];  %%total coordinates
Alist = [GB_atyp'; bulk_atyp'];   %%total list
Clist = [ones(GB_num, 1); zeros(bulk_num,1)]; %chainlist

item = 0;
for type=1:length(numIons)
    s = find(Alist(:)==ORG_STRUC.atomType(type));
    for i=1:length(s)
        T_coor = [T_coor; coord(s(i),:)];
        chanAList = [chanAList; Clist(s(i))];
    end
    T_atyp(1,item+1:item+numIons(type)) = ORG_STRUC.atomType(type);
    item = item + numIons(type);
end
