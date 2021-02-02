function [coor, lat] = Read_SIESTA_Structure()
%0.10033 -3.5268 -1.6877
%3.5629 0.6309 -3.3755
%1.11 -0.31545 1.6877
%6
%1 6 0.33154 0.31666 0.55295
%1 6 0.0076006 0.48422 0.089658
%1 6 0.84599 0.29212 0.83169
%1 6 0.70713 0.57499 0.45019
%1 6 0.85566 0.025895 0.68842
%1 6 0.24248 0.04263 0.52085

coor = [];
lat = zeros(3,3);
if ~exist('OUT.UCELL.ZMATRIX')
   fp = fopen('siesta.STRUCT_OUT');
   count=1;
   
   for i=1:3
       tmp = fgetl(fp);
       lat(i,:) = str2num(tmp);
   end
   tmp = fgetl(fp);
   N_atom = sum(str2num(tmp));
   for i = 1:N_atom
       tmp = str2num(fgetl(fp));
       coor(count,:) = tmp(3:5);
       count = count+1;
   end
   coor = coor - floor(coor);
   fclose(fp);
else %molecule, will read zmatrix+lat
%LatticeConstant 1.0 Ang
%%block LatticeVectors
%          8.300000000       0.000000000       0.000000000
%          0.000000000       8.300000000       0.000000000
%          0.000000000       0.000000000       8.300000000
%%endblock LatticeVectors
%ZM.UnitsLength Ang
%ZM.UnitsAngle rad
%%block Zmatrix
%molecule_cartesian
% 1   0   0   0      0.00000000      0.00000000      0.00000000   0   0   0
% 2   1   0   0      1.10006000      1.49384032      1.24365728   0   1   1
% 2   1   2   0      1.10006000      1.91063300     -0.55870990   0   0   1
% 2   1   2   3      1.10006000      1.91063300      2.09439500   0   0   0
% 2   1   2   3      1.10006000      1.91063300      4.18879000   0   0   0

   fp = fopen('OUT.UCELL.ZMATRIX');
   running = 1;
   count = 1;
   dolat = 0;
   docoor = 0;
   while running
         a = fgetl(fp);
      if findstr(a, '%block LatticeVectors')
          dolat = 1;
      if findstr(a, '%endblock LatticeVectors')
          dolat=0;
      elseif findstr(a, '%block Zmatrix')
          docoor = 1;
      elseif findstr(a, '%endblock Zmatrix')
          running=0;
          docoor = 0;
      elseif ~findstr(a, 'molecule_cartesian')
          if dolat==1
             optlat(1,:) = str2num(a);
             optlat(2,:) = str2num(fgetl(fp));
             optlat(3,:) = str2num(fgetl(fp));
          elseif  (docoor==1)  
             tmp = str2num(a);
             coor(count,:) = tmp(5:7);
             count = count + 1;
          end
      end
   end
   fclose(fp);
end
end
