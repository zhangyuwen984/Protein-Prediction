function lat_OK = latticeCheck(Lattice)

% USPEX Version 6.2
% Change: added 3rd constrain
% Returns 1 in case the lattice is OK, 0 - otherwise.

% this function checks whether the given lattice fulfills the hard
% constraints for lattices. There are three sorts of hard constraints.
% 1) A minimal angle between any two vectors defining the lattice.
% 2) A minimal distance between any two planes; was lately changed to simply lattice vector length
% 3) A minimal angle between the vector defining the lattice and the
% diagonal of the parallelogram formed by other 2 vectors defining the lattice


global ORG_STRUC


% lat_OK
lat_OK = 1;

% in this function both the matrice and the parameter represatation is
% required, so first we prepare the two
if length(Lattice)==3
    angLattice = latConverter(Lattice);
else
    angLattice = Lattice;
    Lattice = latConverter(Lattice);
end

% calculate the angles is degrees
angles = angLattice(4:6)*180/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is constraint 1)
% check whether none of the angles is two small. Note that the problem is
% symmetric around 90 degrees, that's the reason for 180-minAngle.

if (isempty(find(angles<ORG_STRUC.minAngle)) + isempty(find(angles>(180-ORG_STRUC.minAngle))))~=2
    lat_OK = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is constraint 2)
% We receive the distance by dividing the volume (found by det(lattice)) by the area of any two vectors

vol = abs(det(Lattice));
dist = zeros(3,1);
dist(1) = vol/(angLattice(1)*angLattice(2)*sin(angLattice(6)));
dist(2) = vol/(angLattice(1)*angLattice(3)*sin(angLattice(5)));
dist(3) = vol/(angLattice(2)*angLattice(3)*sin(angLattice(4)));

%if ~isempty(find(dist<minVectorLength))
if ~isempty(find(angLattice(1:3)<ORG_STRUC.minVectorLength))
    lat_OK = 0;
end

if ~isreal(Lattice)
    lat_OK = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is constraint 3)
% check whether none of the angles between the vector defining the lattice and the
% diagonal of the parallelogram formed by other 2 vectors defining the
% lattice is two small.

if lat_OK % if it's not 0 it means there are no 0-length vectors
    a_bc = Lattice(1,1)*(Lattice(2,1)+Lattice(3,1))+Lattice(1,2)*(Lattice(2,2)+Lattice(3,2))+Lattice(1,3)*(Lattice(2,3)+Lattice(3,3));
    ab_c = Lattice(3,1)*(Lattice(2,1)+Lattice(1,1))+Lattice(3,2)*(Lattice(2,2)+Lattice(1,2))+Lattice(3,3)*(Lattice(2,3)+Lattice(1,3));
    b_ca = Lattice(2,1)*(Lattice(1,1)+Lattice(3,1))+Lattice(2,2)*(Lattice(1,2)+Lattice(3,2))+Lattice(2,3)*(Lattice(1,3)+Lattice(3,3));
    
    % |b+c|^2=|b|^2+|c|^2+2|b||c|cos(bc)
    anglesDiag1(1) = acos(a_bc/(angLattice(1)*sqrt(angLattice(2)^2+angLattice(3)^2+2*angLattice(2)*angLattice(3)*cos(angLattice(4)))));
    anglesDiag1(2) = acos(b_ca/(angLattice(2)*sqrt(angLattice(1)^2+angLattice(3)^2+2*angLattice(1)*angLattice(3)*cos(angLattice(5)))));
    anglesDiag1(3) = acos(ab_c/(angLattice(3)*sqrt(angLattice(2)^2+angLattice(1)^2+2*angLattice(2)*angLattice(1)*cos(angLattice(6)))));
    anglesDiag = anglesDiag1(1:3)*180/pi;
    
    if (isempty(find(anglesDiag<ORG_STRUC.minDiagAngle)) + isempty(find(anglesDiag>(180-ORG_STRUC.minDiagAngle))))~=2
        lat_OK = 0;
    end
end
if ORG_STRUC.dimension==2
    lat_OK = 1
end
if ORG_STRUC.constLattice
    lat_OK = 1;
end
