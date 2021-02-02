function [lat_n, coor_n, order_n, numIons_n]= supercell_301(lat, coor, order, numIons, cell1, cell2)

%%This program is to covert the cell1 to the cell2
%%The basic idea is like this, let take the example of conversion from 4 to 5
%%firstly create super cell of 20, and then cut it to 5, and take one region randomly.
       lat_p = latConverter(lat);
       [b,I] = sort(lat_p(1:3));
       coor_n =[];
       order_n =[];
       ntype = size(numIons,2);
       numIons_n = zeros(1, ntype);
       lat_n = lat;
       lat_n(I(1),:)=lat(I(1),:)*cell2/cell1;
%%%%%1:create super cell along the shortest direction
       supercell = lcm(cell1,cell2);
       num1=supercell/cell1;
       count1=0;
       count2=0;

       for i=1:ntype
           for j = 1:numIons(i)
                count1=count1+1;
                for k=0:num1-1
                    count2=count2+1;
                     coor_s(count2,:)= coor(count1,:);
                     coor_s(count2,I(1))= coor(count1,I(1))/num1 +k/num1;
                    order_s(count2)=order(count1);       
                     type_s(count2) = i;
                end
           end          
       end
%%%%%%2: chose the region
      num2=supercell/cell2;
      coor_s=coor_s-floor(coor_s);
%%%%%%%%%%%%%%how many atoms in each region
      pick = RandInt(1,1,[1,num2]);
      count3 = 1;
      for i=1:count2
	  lower = (pick-1)/num2;
	  upper = pick/num2;
	  tmp = coor_s(i,(I(1)));
          if (tmp > lower) & (tmp < upper)
	      coor_n(count3,:) = coor_s(i,:);
              coor_n(count3,I(1))=(tmp-lower)/(upper-lower);	
	      order_n(count3) = order_s(i);
	      numIons_n(type_s(i)) = numIons_n(type_s(i)) + 1;
              count3=count3+1;			  
           end
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
