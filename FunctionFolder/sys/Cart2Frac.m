function frac = Cart2Frac(cart, lat)
%X, Y, and Z are converted into fractional crystallographic coordinates (x,y,z) 
%in order to perform crystallographic operations.
%Inversely, geometric computations are more easily performed in Cartesian space. 
%In orthonormal systems (cubic, tetragonal, and orthorhombic), the coordinate 
%transformation reduces to a simple division of the coordinate values by cell constants. 
%In the case of a generic oblique crystallographic system, 
%the transformation is described by a matrix operation:
%http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm
%Created by Qiang Zhu (2014/01/02)

a     = lat(1);
b     = lat(2);
c     = lat(3);
alpha = lat(4);
beta  = lat(5);
gamma = lat(6);
v = a*b*c*(1-cos(alpha)^2 - cos(beta)^2 -cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma))^(1/2);
%M = [a b*cos(gamma) c*cos(beta); ...
%     0 b*sin(gamma) c*(cos(alpha)-cos(beta)*cos(gamma)/sin(gamma));...
%     0 0            c*v/sin(gamma)];
M = [ 1/a, -cos(gamma)/(a*sin(gamma)), (b*cos(gamma)*c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)-b*c*cos(beta)*sin(gamma))/v;...
      0,    1/(b*sin(gamma)),           -a*c*(cos(alpha)-cos(beta)*cos(gamma))/(v*sin(gamma)); ...
      0,            0,                  a*b*sin(gamma)/v];
frac = (M*cart')';

