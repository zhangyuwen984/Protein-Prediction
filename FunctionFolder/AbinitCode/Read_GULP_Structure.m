function [coor, lat] = Read_GULP_Structure()
%This routine is to read crystal structure from GULP output
%File: optimized.structure
%fractional for bulk
%cartesian for surface
%Last updated by Qiang Zhu (2013/10/03)

fp = fopen('optimized.structure');
running = 1;
count = 1;
dolat = 0;
docoor = 0;
while running
    a = fgetl(fp);
    if findstr(a, 'cell')
           dolat = 1;
%make sure it works with octave
	elseif ~(isempty(findstr(a, 'fractional')) & isempty(findstr(a, 'cartesian')))
           docoor = 1;
    elseif dolat==1
           optlat = str2num(a);
           if length(optlat)>3
              optlat(4:6) = optlat(4:6)*pi/180;
              lat = latConverter(optlat);
           else
              lat = 0; %surface
           end
           dolat=0;
     elseif  docoor == 1
           if findstr(a, 'core')
              tmp = str2num(a(12:40));
              coor(count,:) = tmp;
              count = count + 1;
           else
               running = 0;
               docoor=0;
           end
     end

end
fclose(fp);
coor = coor - floor(coor);
