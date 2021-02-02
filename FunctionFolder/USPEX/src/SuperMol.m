function [Molecules1, numMols1, lat1] = SuperMol(Molecules, numMols, lat, splitInto)

total = prod(splitInto);
numMols1=numMols*total;

lat1=lat;
for i=1:3
    lat1(i,:)=lat(i,:)*splitInto(i);
end


item = 1;
for i=1:splitInto(1)
   for j=1:splitInto(2)
       for k=1:splitInto(3)
           ABS(item,:)=lat(1,:)*(i-1)+lat(2,:)*(j-1)+lat(3,:)*(k-1);
           item = item + 1;
       end
   end
end

item = 0;
for i=1:sum(numMols)
    for j = 1:total
        Molecules1(item+j) = Molecules(i);
        Molecules1(item+j).ZMATRIX(1,:) = Molecules(i).ZMATRIX(1,:) + ABS(j,:);
        for k = 1:size(Molecules(i).MOLCOORS,1)
            Molecules1(item+j).MOLCOORS(k,:) = Molecules(i).MOLCOORS(k,:) + ABS(j,:);
        end
    end
    item = item + total;
end
