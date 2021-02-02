function CellList=findcell(N)
CellList=[];
item = 1;
for i=1:N
   for j=1:N
       if i*j<=N
          CellList(item,:)=[i,j];
          item = item + 1;
       end
   end
end
