function [coor, lat] = Read_DMACRYS_Structure()
%to Read GULP

fp = fopen('fort.13');
running = 1;
count = 1;

while running
    a = fgetl(fp);
    if findstr(a, 'CELL')
           optlat = str2num(a(6:76));
           optlat(4:6) = optlat(4:6)*pi/180;
           lat = latConverter(optlat);
    elseif findstr(a, 'ATOM') 
           tmp = str2num(a(20:58));
           coor(count,:) = tmp;
           count = count + 1;
    elseif findstr(a, 'END') 
           running = 0;
    end
end
fclose(fp);
