function [Zmatrix, MOLCOORS] = RotInertia(coords, format, num_optFlags, molRot)

global ORG_STRUC

STDMOL = ORG_STRUC.STDMOL;
MOLCOORS = coords;
[a,b] = PrincipleAxis(coords);  %find the principle roation axis

%b(1,1) =0 for linear chain molecules, thus two rotational variables
if abs(b(1,1)) < 0.0001 
   Ref = b(2,2);
else
   Ref = b(1,1);
end

for loop = 1:3
    if b(loop, loop) > 0.0001
       angle = (rand-0.5)*pi/2*Ref/b(loop,loop);
       MOLCOORS = Rotate_rigid_body(mean(MOLCOORS), a(:,loop)', MOLCOORS, angle);
       %disp(['Rot Angle: ' num2str(loop) '-' num2str(angle)])
    end
end

MOLCOORS = bsxfun(@plus, MOLCOORS, rand(1,3)-0.5); %Random Translation
Zmatrix=NEW_coord2Zmatrix(MOLCOORS, format);

if num_optFlags > 0
    goodRot = 0;
    for i = 1:size(coords, 1)
        radiu(i) = str2num(covalentRadius(ORG_STRUC.atomType(STDMOL(molRot).types(i))));
    end

    while ~goodRot
        for ind = 1: num_optFlags
            i1 = STDMOL(molRot).flex_dihedral(ind,1);
            i2 = STDMOL(molRot).flex_dihedral(ind,2);
            Zmatrix(i1,i2) = Zmatrix(i1,i2) + (rand-0.5)*pi;
        end
        MOLCOORS = NEW_ZMATRIXCOORD(Zmatrix, format);
        Zmatrix  = NEW_coord2Zmatrix(MOLCOORS, format);
        CN = find_pair(MOLCOORS, radiu);
        if isequal(CN, STDMOL(molRot).CN)
           goodRot = 1;
        end
    end
end

