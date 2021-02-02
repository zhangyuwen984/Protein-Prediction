function GB_order=reduce_order(order, chanAList)

     GB_order = [];

     count = 0;
     for i = 1:length(chanAList)
         if chanAList(i)==1
            count = count + 1;
            GB_order(count)= order(i);
         end
     end
