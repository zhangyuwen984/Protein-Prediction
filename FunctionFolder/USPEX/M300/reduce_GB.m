function GB_coor=reduce_GB(lattice,coordinates,chanAList, thickness, GB_numIons)


     lat=latConverter(lattice);
     
     bulk_thickness=lat(3)-thickness;

     GB_coor = zeros(sum(GB_numIons),3);

     count = 0;
     for i = 1:length(chanAList)
         if chanAList(i)==1
            count = count + 1;
            GB_coor(count,:)= coordinates(i,:);
            GB_coor(count,3)= GB_coor(count,3)*lat(3)-bulk_thickness/2;
            GB_coor(count,3)= GB_coor(count,3)/(lat(3)-bulk_thickness);
          %   if GB_coor(count,3) < 2 % go to upper surface
          %      GB_coor(count,3) = GB_coor(count,3)/thickness +1 ;
          %   else
          %      GB_coor(count,3) = (GB_coor(count,3)- bulk_thickness)/thickness ;
          %   end

            %we allow the GB atoms go to +/- 1A of substrate
          %  if (GB_coor(count, 3) < 0) & GB_coor(count,3)*thickness < -1
          %        GB_coor(count, 3) = GB_coor(count, 3) + 1;
          %        disp('too small, reset')
          %  elseif (GB_coor(count, 3) > 1) & ((GB_coor(count,3)-1)*thickness > 1)
          %        GB_coor(count, 3) = GB_coor(count, 3) - 1;
          %        disp('too big, reset')
          %  end
         end
     end
    
     
