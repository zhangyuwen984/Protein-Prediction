function newCoords = Add_Atom(lat, coor, type, NAdd, atomType)

%Delete atoms according to the height in c-axis
if sum(NAdd) == 0
    newCoords = coor;
else
    newCoords = [];
    adatom_type = [];
    adatom_coord = [];
    Occupied=zeros(1,length(atomType));
    
    for i=1:sum(NAdd)
        goodpick = 0;
        while ~goodpick
    %%%%try to pick the atom to be added randomly
          t=RandInt(1,1,[1,length(atomType)]);
          if Occupied(t) < NAdd(t)+0.5
             goodpick=1;
             [x,y]=choose_docking([adatom_coord;coor],lat,[adatom_type;type']',atomType(t));
             zmin = 1/lat(3,3);
             z = surface_adatom_docking(lat,zmin,[adatom_coord;coor],[adatom_type;type'],[x,y],atomType(t));
             adatom_coord=[adatom_coord;[x,y,z]];
             adatom_type = [adatom_type;atomType(t)];
             Occupied(t)=Occupied(t)+1;
          end
        end
    end
    adatom_type = [adatom_type;type'];
    adatom_coord = [adatom_coord;coor];
    for ind = 1: length(atomType)
        s = find(adatom_type == atomType(ind));
        newCoords = cat(1,newCoords, adatom_coord(s,:));
    end
end
