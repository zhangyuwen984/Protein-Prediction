function [coor, lat] = Read_CASTEP_Structure()

coor = [];

handle = fopen('cstp-out.cell');
while 1
    lat_tmp = fgetl(handle);
    if lat_tmp == -1
        break
    end
    if ~isempty(findstr(lat_tmp, 'ANG'))
	for LL = 1 : 3
	    lat_tmp = fgetl(handle);
	    lat_1 = sscanf(lat_tmp,'%g %g %g');
	    lat(LL,1) = lat_1(1);
	    lat(LL,2) = lat_1(2);
	    lat(LL,3) = lat_1(3);
	end
	% lat is 3x3 matrix now.
    end
end
status = fclose(handle);

handle = fopen('cstp-out.cell');
stop = 0;
while 1
    tmp = fgetl(handle);
    row_num = 0;
    if (~isempty(findstr(tmp, 'positions_frac')))
	while 1
	    row_num   = row_num + 1;
	    tmp       = fgetl(handle);
	    row_array = sscanf(tmp, '%*s %g %g %g');
	    if (~isempty(findstr(tmp, 'positions_frac')))
		stop = 1;
		break
	    else
		coor(row_num, 1) = row_array(1);
		coor(row_num, 2) = row_array(2);
		coor(row_num, 3) = row_array(3);
	    end
	end
    end
    if tmp == -1 | stop == 1
	break
    end
end
status = fclose(handle);

coor = coor - floor(coor);
