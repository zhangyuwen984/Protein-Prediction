function [dist]=cosineDistance(matrA, matrB, weight)

% we assume that both matrices have the same dimensions

% dist = 1-(dot(vecA, vecB))/(sqrt(dot(vecA, vecA)*dot(vecB,vecB)));

%[r, c] = size(matrA);
%
%coef1 = 0;
%coef2 = 0;
%coef3 = 0;
%
%for i = 1:r
%  for j = 1:c
%   coef1 = coef1 + matrA(i,j)*matrB(i,j)*weight(i);
%   coef2 = coef2 + matrA(i,j)*matrA(i,j)*weight(i);
%   coef3 = coef3 + matrB(i,j)*matrB(i,j)*weight(i);
%  end
%end

%let's write it in a concise way
weight = diag(weight);
coef1 = sum(sum(weight*( matrA.*matrB )));
coef2 = sum(sum(weight*( matrA.*matrA )));
coef3 = sum(sum(weight*( matrB.*matrB )));

dist = (1-coef1/(sqrt(coef2*coef3)))/2;
if abs(dist-1) < 0.000001
   dist = 0.99999;
   disp(['Fingerprint might be wrong !!!']);
end

