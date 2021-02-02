function [coor, composition] = getSurface(coordinate, numIons, lat)

global ORG_STRUC

item=1;
k=1;
composition=zeros(size(numIons));
coor=[];
for i=1:length(numIons)
    for j=1:numIons(i)
        if coordinate(k,3)*lat(3,3) > ORG_STRUC.bulk_lat - ORG_STRUC.thicknessB
           coor(item,:)=coordinate(k,:);
           item=item+1;
           composition(i)=composition(i)+1;
        end
    k=k+1; %coordinates
    end
end
