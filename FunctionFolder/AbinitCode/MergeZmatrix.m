function [MOLCOOR, height] = MergeZmatrix(Zmatrix, format, bound, type)
%This utility is to merge two monomers. The basic idea is to connect two boundary atoms together.
% M11+++++M12 --- M21+++++M22
  N_atom   = length(type);
  Zmatrix1 = Zmatrix;
  coor1    = NEW_ZMATRIXCOORD(Zmatrix,format);

%all we need to do is to redefine the first 3 components:
line = find(bound(:)==1);
if length(line) ~= 2
   disp(['we only support 2 connecting atoms in 110']);
   quit
end

% M11+++++M12 --- M21+++++M22
   connect = line(2);

% to calculate the coordination number of M12?
count=1;
MATRIX(count,:)=coor1(connect,:);
for i=1:N_atom
    if i~=connect
          dist=norm(coor1(i,:)-coor1(connect,:));
          radii1 = str2num(covalentRadius(type(i)));
          radii2 = str2num(covalentRadius(type(connect)));
       if dist<1.2*(radii1+radii2)
          count = count + 1;
          MATRIX(count,:) = coor1(i,:);
       end
    end
end

%to locate M21
    Zmatrix2 = Zmatrix1;
if count == 4 %4-coordination
    vector = (MATRIX(2,:) + MATRIX(3,:) + MATRIX(4,:))/3 -MATRIX(1,:);
    vector = vector/norm(vector);
    Zmatrix2(1,:)= MATRIX(1,:)-1.5*vector;
elseif count == 3 %3-coordination
    Zmatrix2(1,:)=3*MATRIX(1,:)-MATRIX(2,:)-MATRIX(3,:);
elseif count == 2 % 2-coordination
    Zmatrix2(1,:)=2*MATRIX(1,:)-MATRIX(2,:);
end
   coor2 = NEW_ZMATRIXCOORD(Zmatrix2, format);

%to rotate around the center, we do it twice to make sure it is well done.
    center1 = mean(coor1);
    center2 = mean(coor2);
    center = 0.5*(center1+center2);
    coor1 = bsxfun(@minus, coor1, center);  %Deviations from Mean
    vector1 = center1-center;
    vector1 = vector1/norm(vector1);
    R_axis = cross(vector1,[0 0 1]);
    R_axis = R_axis/norm(R_axis);
    angle  = -acos(dot([0 0 1], vector1));
    MOLCOOR  = Rotate_rigid_body([0 0 0], R_axis, coor1, angle);
    height = norm(center1-center2);
