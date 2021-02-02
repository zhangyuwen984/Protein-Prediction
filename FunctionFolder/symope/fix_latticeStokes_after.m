function [lat, coord] = fix_latticeStokes_after(nsym, oldLat, oldCoord)
% fix some Stokes 'problems', for example non conventional choice of primitive monoclinic lattice

lat = oldLat;
coord = oldCoord;
sGroup = spaceGroups(nsym); % space group's standard symbol

if (nsym > 2) & (nsym < 16) & (sGroup(1) == 'P') % monoclinic, primitive
    % 8 7 10 90 90 100 => 7 10 8 90 100 90
    lat(5) = oldLat(6); % non-90 angle is gamma for Stokes and beta conventionally
    lat(6) = 90;
    lat(4) = 90;
    lat(1) = oldLat(2); % have to swap axes as well
    lat(2) = oldLat(3); % have to swap axes as well
    lat(3) = oldLat(1); % have to swap axes as well
    coord(:,1) = oldCoord(:,2); % have to swap coordinates as well
    coord(:,2) = oldCoord(:,3);
    coord(:,3) = oldCoord(:,1);
end
