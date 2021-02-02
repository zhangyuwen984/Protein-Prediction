function [coor, lat] = Read_LAMMPS_Structure()
fp=fopen('out.xyz','r');
while 1
    tline = fgetl(fp);
    if ~ischar(tline)
    fclose(fp);
    break;
    end
    tmp=sscanf(tline,'%d');
    NCoords=tmp;
    tline = fgetl(fp);
        for i=1:NCoords
            tline = fgetl(fp);
            tmp=sscanf(tline,'%f %f %f %f');
            coords(i,:)=tmp(2:4);

        end

end
fp=fopen('lammps.out','r');
while 1
    tline = fgetl(fp);
    if ~isempty(findstr(tline, 'Lattice parameters'))
        for i=1:3
            tline=fgetl(fp);
            cell(i,:)=sscanf(tline,'%f %f %f');
        end
    end
    if ~ischar(tline)
    fclose(fp);
    break;
    end
end
coords=coords*(cell^-1);
lat=cell;
coor=coords;
end
