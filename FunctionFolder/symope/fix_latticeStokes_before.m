function [lat] = fix_latticeStokes_before(nsym, oldLat)
% fix some Stokes 'problems', for example non conventional choice of primitive monoclinic lattice

lat = oldLat;
sGroup = spaceGroups(nsym); % space group's standard symbol

if (nsym > 2) & (nsym < 16) & (sGroup(1) == 'P') % monoclinic, primitive
    % 7 10 8 90 100 90 => 8 7 10 90 90 100
    %lat(6) = oldLat(5); % non-90 angle is gamma for Stokes and beta conventionally
    %lat(5) = 90;
    %lat(4) = 90;
    %lat(3) = oldLat(2); % have to swap axes as well
    %lat(1) = oldLat(3); % have to swap axes as well
    %lat(2) = oldLat(1); % have to swap axes as well
    lat(1) = oldLat(3); % have to swap axes as well
    lat(3) = oldLat(1); % have to swap axes as well
end
