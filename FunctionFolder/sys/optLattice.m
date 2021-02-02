function [optcoor,optlat] = optLattice(coor,lattice)

% idea is as follow - if any lattice vector has projection onto any other lattice vector
% that is greater than half of length of this vector, we can reoptimise the shape
% i. e. if |a*b|/|b| > |b|/2 then a_new = a - ceil(|a*b|/|b|^2)*sign(a*b)*b

v1 = lattice(1,:);
v2 = lattice(2,:);
v3 = lattice(3,:);

%coor = coor*lattice; % now coordinates are absolute

flag = 1;
step = 0;
while flag
    flag = 0;
    [v1, flag] = reoptimizeVector(v1,v2,flag);
    [v1, flag] = reoptimizeVector(v1,v3,flag);
    [v2, flag] = reoptimizeVector(v2,v1,flag);
    [v2, flag] = reoptimizeVector(v2,v3,flag);
    [v3, flag] = reoptimizeVector(v3,v1,flag);
    [v3, flag] = reoptimizeVector(v3,v2,flag);
    
    if step > 100
        break;    % so that the cycle isn't infinite
    end
    
    [v1, flag] = reoptimizeVector(v1,v2+v3,flag);
    [v2, flag] = reoptimizeVector(v2,v1+v3,flag);
    [v3, flag] = reoptimizeVector(v3,v1+v2,flag);
    
    if flag
        step = step + 1;
    end
end

optlat(1,:) = v1;
optlat(2,:) = v2;
optlat(3,:) = v3;

if (det(optlat) < 0.0000001) & (sum(coor(:)) == 0)
    tempcoor = coor;
    optlat = [1 0 0; 0 1 0; 0 0 1];
else
    tempcoor = coor/optlat;
end

optcoor = tempcoor - floor(tempcoor); % now coordinates are fractional
optcoor = optcoor*optlat;             % now absolute again for molecules

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, flag] = reoptimizeVector(v1,v2,flag1)

flag = flag1;
v = v1;
if abs(dot(v1,v2)) > (norm(v2)/2)
    try
        v1_trial = v1 - ceil(abs(dot(v1,v2))/(norm(v2)^2))*sign(dot(v1,v2))*v2;
    catch
        v1_trial = v1;
    end
    if norm(v1_trial) < norm(v1)
        v = v1_trial;
        flag = 1;
    end
end
