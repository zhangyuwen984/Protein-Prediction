function Zmatrix= NEW_coord2Zmatrix(coords, format)

pi=3.1415926;
if(size(coords,1)==1)
 Zmatrix=coords;
else 

  Zmatrix = coords(1,:);
  coords(:,1)=coords(:,1)-coords(1,1);
  coords(:,2)=coords(:,2)-coords(1,2);
  coords(:,3)=coords(:,3)-coords(1,3);
  Zmatrix(2,1)=(sum((coords(2,:)-coords(1,:)).^2))^(0.5);

   if coords(2,3)==0
      Zmatrix(2,2)=pi/2;
   else
       Zmatrix(2,2)=acos(coords(2,3)/Zmatrix(2,1));
   end
   if coords(2,2)==0
       Zmatrix(2,3)=0;
   else
       Zmatrix(2,3)=atan2(coords(2,2),coords(2,1));
   end
for ind = 3:size(coords,1)
   Zmatrix(ind,1)=(sum ((coords(ind,:)-coords(format(ind,1),:)).^2 ))^(0.5);
   vectorsCrutch2= (coords(format(ind,1),:)-coords(format(ind,2),:));
   vectorsCrutch1=(coords(ind,:)-coords(format(ind,1),:));
   if ind==3
       vectorsCrutch3=[0,0,-1];
   else
   vectorsCrutch3= coords(format(ind,2),:)-coords(format(ind,3),:);
   end
    Zmatrix(ind,2)= acos(dot(vectorsCrutch1,-vectorsCrutch2)/((sum(vectorsCrutch1.^2)).^(0.5) * sum(vectorsCrutch2.^2).^(0.5)));
    normalVEC1 = cross(vectorsCrutch2,vectorsCrutch3);
    normalVEC2 = cross(vectorsCrutch1,vectorsCrutch2);
    torsang = atan2(norm(vectorsCrutch2)*dot(vectorsCrutch1,normalVEC1),dot(normalVEC2,normalVEC1));
    Zmatrix(ind,3)=-torsang;
    if abs(torsang+pi)<0.0001
        Zmatrix(ind,3)=pi;
    end
end

end
