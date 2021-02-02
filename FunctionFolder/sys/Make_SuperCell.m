function [S_coor, S_lat, N_Size] = Make_SuperCell(coor, lat, Dim, flag)

%To make the supercell
%if Dim is a number (i), we make i*i*i
%if Dim is a 1*2 vector ([i,j]), we make i*j*1
%if Dim is a 1*3 vector ([i,j,k]), we make i*j*k
%flag: 1, Write the replication atom by atom 
%flag: 0, Write the replication cell by cell 
Size = [1; 1; 1];
if length(Dim) == 1
    Size= [Dim; Dim; Dim];
elseif length(Dim) == 2
    Size= [Dim(1); Dim(2); 1];
elseif length(Dim) == 3
    Size= [Dim(1); Dim(2); Dim(3)];
else
    disp(['INPUT Dimension ' num2str(Dim) ' is incompatible']);
    disp(['will return 1*1*1 cell']);
end

N_Size = prod(Size);
N_atom = size(coor,1);

if N_Size == 1
    S_coor = coor;
    S_lat  = lat;
else
    Matrix = SuperMatrix(0, Size(1)-1, 0, Size(2)-1, 0, Size(3)-1);
    %X1 = [0:Size(1)-1];
    %Y1 = [0:Size(2)-1];
    %Z1 = [0:Size(3)-1];
    %[X2,Y2,Z2] = meshgrid(X1,Y1,Z1);
    %X3=reshape(X2,1,[]); 
    %Y3=reshape(Y2,1,[]); 
    %Z3=reshape(Z2,1,[]);
    %Matrix = [X3; Y3; Z3]';

    if flag == 0 %cell by cell
       S_coor   = repmat(coor, [N_Size,1]);
       %  ATOM: a1 b1 c1   ---->    a1 b1 c1
       %        a2 b2 c2   ---->    a2 b2 c2
       %                            a1 b1 c1
       %                            a2 b2 c2
       tmp      = repmat(Matrix, [1, N_atom]);
       %  Trans: 0 0 0           0 0 0 0 0 0
       % Matrix: 0 0 1  ------>  0 0 1 0 0 1
       %         0 1 0           0 1 0 0 1 0

       S_Matrix = reshape(tmp',[3,N_atom*N_Size])';
       %  Trans: 0 0 0           0 0 0 
       % Matrix: 0 0 1  ------>  0 0 0 
       %         0 1 0           0 0 1 
       %         0 0 0           0 0 1 
       %         0 0 1           0 1 0 
       %         0 1 0           0 1 0 

    else %atom by atom, the other way
       tmp      = repmat(coor, [1, N_Size]);
       S_coor   = reshape(tmp',[3,N_atom*N_Size])';
       S_Matrix = repmat(Matrix, [N_atom,1]);
    end

    S_coor  = S_coor + S_Matrix;
    S_coor  = S_coor * lat;
    S_lat   = repmat(Size,[1,3]).*lat;
    S_coor  = S_coor / S_lat;
end
