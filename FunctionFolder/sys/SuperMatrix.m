function Matrix = SuperMatrix(xmin, xmax, ymin, ymax, zmin, zmax)

%This is a small utlity to quickly generate a 3d matrix series 
%such as [1 0 0; 0 1 0; 0 0 1; ........]
%The functional is usually used to create super cell
%Example
%INPUT:  [0 1 0 1 0 1]
%OUTPUT: [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]

%X1 = [xmin:xmax];
%Y1 = [ymin:ymax];
%Z1 = [zmin:zmax];
%[X2,Y2,Z2] = meshgrid(X1,Y1,Z1);
%X3=reshape(X2,1,[]);
%Y3=reshape(Y2,1,[]);
%Z3=reshape(Z2,1,[]);
%Matrix = [X3; Y3; Z3]';
%ToDelete = all(Matrix==0, 2);
%Matrix(ToDelete,:) =[]; %[0 0 0]

%-------------Equivalent to the following
Matrix = [];
for i = xmin:1:xmax
    for j = ymin:1:ymax
        for k = zmin:1:zmax
            Matrix = [Matrix; [i,j,k]];
        end
    end
end
