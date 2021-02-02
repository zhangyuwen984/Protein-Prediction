function [a,b] = PrincipleAxis(AbsoluteCoord)

AbsoluteCoord = bsxfun(@minus, AbsoluteCoord, mean(AbsoluteCoord)); % Matrix-Vector;

Inertia = zeros(3,3); % moment of inertia tensor
Inertia(1,1) = sum(AbsoluteCoord(:,2).^2 + AbsoluteCoord(:,3).^2);
Inertia(2,2) = sum(AbsoluteCoord(:,1).^2 + AbsoluteCoord(:,3).^2);
Inertia(3,3) = sum(AbsoluteCoord(:,1).^2 + AbsoluteCoord(:,2).^2);
Inertia(1,2) = -sum(AbsoluteCoord(:,1).* AbsoluteCoord(:,2));
Inertia(2,3) = -sum(AbsoluteCoord(:,2).* AbsoluteCoord(:,3));
Inertia(3,1) = -sum(AbsoluteCoord(:,3).* AbsoluteCoord(:,1));
Inertia(2,1) = -sum(AbsoluteCoord(:,1).* AbsoluteCoord(:,2));
Inertia(3,2) = -sum(AbsoluteCoord(:,2).* AbsoluteCoord(:,3));
Inertia(1,3) = -sum(AbsoluteCoord(:,3).* AbsoluteCoord(:,1));
%Vector*Vector
[a,b] = eig(Inertia);   % maximal moment of inertia - for Z axis, minimum - for X

