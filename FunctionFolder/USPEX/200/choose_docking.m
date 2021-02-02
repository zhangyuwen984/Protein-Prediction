function [x,y]=choose_docking(candidate,lat, typelist, type)

%%%%%%%%%This routine is used to compensate the atom if some add atoms are removed from the network during relaxation%%%%%%%%%%%
%%%%%%%%pick the atoms at top layer
candidate1 = candidate*lat;
lat1(1,1)=lat(1,1);
lat1(1,2)=lat(1,2);
lat1(2,1)=lat(2,1);
lat1(2,2)=lat(2,2);
lat = latConverter(lat);
item = 0;
coor=[];
for i = 1:size(candidate,1)
    if candidate1(i,3) > max(candidate1(:,3))-1
        item = item + 1;
        coor(item,:) = candidate(i,:);
        types(item) = typelist(i);
    end
end
%%%%%%%%%start to mesh%%%%%%%%%%%
maxdist = 0;
mesh1 = floor(lat(1))*1.25;
mesh2 = floor(lat(2))*1.25;
count = 0;
goodmesh = [];
bestmesh = [];
for i = 1:mesh1
    for j = 1:mesh2
        sumdist = 0;
        good = 1;
        for k= 1:item
           dist = distance(lat1, [i/mesh1,j/mesh2], coor(k,1:2) );
           radii1 = str2num(covalentRadius(ceil(type)));
           radii2 = str2num(covalentRadius(ceil(types(k))));
           sumdist = sumdist + dist;
           if dist < (radii1+radii2)
               good = 0;
           end
        end

        if sumdist > maxdist
           maxdist = sumdist;
           bestmesh = [i/mesh1,j/mesh2];
        end

        if good
           count = count + 1;
           goodmesh(count,:)=[i/mesh1,j/mesh2];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%randomly pick 
if count < 0.5
   x = bestmesh(1);
   y = bestmesh(2);
else
   num = RandInt(1,1,[1,count]);
   x = goodmesh(num,1);
   y = goodmesh(num,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function dist=distance(lattice, coord1, coord2)
     check = coord1-coord2;
     for index=1:2
        check(index)=check(index)-round(check(index));
     end
        dist = sqrt(sum((check*lattice).^2));
