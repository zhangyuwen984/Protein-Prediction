function output = Rotate_rigid_body(point, R_axis, input, angle)
%To Rotate some angle along a given axis, relative to a reference point
r = bsxfun(@minus, input, point);

mu = R_axis/norm(R_axis);  %normalized vector
%http://en.wikipedia.org/wiki/Rotation_matrix: Rotation matrix from axis and angle
M(1,:)=[cos(angle)+mu(1)^2*(1-cos(angle)), mu(1)*mu(2)*(1-cos(angle))-mu(3)*sin(angle), mu(1)*mu(3)*(1-cos(angle))+mu(2)*sin(angle)];
M(2,:)=[mu(1)*mu(2)*(1-cos(angle))+mu(3)*sin(angle), cos(angle)+mu(2)^2*(1-cos(angle)), mu(2)*mu(3)*(1-cos(angle))-mu(1)*sin(angle)];
M(3,:)=[mu(3)*mu(1)*(1-cos(angle))-mu(2)*sin(angle), mu(3)*mu(2)*(1-cos(angle))+mu(1)*sin(angle), cos(angle)+mu(3)^2*(1-cos(angle))];
output = bsxfun(@plus, r*M, point);

