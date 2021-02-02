function  Write_LAMMPS(Ind_No)
% created by Xiao Dong

global POP_STRUC
global ORG_STRUC

Step = POP_STRUC.POPULATION(Ind_No).Step;

numIons= POP_STRUC.POPULATION(Ind_No).numIons;
try
[nothing, nothing] = unix(['cat lammps.in_' num2str(Step) ' > lammps.in']);
catch
   error = ['param file is not present for step' num2str(Step)];
   save ([ ORG_STRUC.resFolder '/ERROR_param.txt'],'error')
   quit
 end
fp = fopen('coo', 'w');
fprintf(fp, 'USPEX\n');
fprintf(fp,'\n');
fprintf(fp,'%d atoms\n',sum(numIons));
fprintf(fp,'%d atom types\n',length(numIons));
fprintf(fp,'\n');
lat1 = latConverter(POP_STRUC.POPULATION(Ind_No).LATTICE);
lat = latConverter(lat1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modified the lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coor=POP_STRUC.POPULATION(Ind_No).COORDINATES;
pos=coor*lat;

lat(2,1)=lat(2,1)-round(lat(2,1)/lat(1,1))*lat(1,1);
lat(3,:)=lat(3,:)-round(lat(3,2)/lat(2,2))*lat(2,:);
lat(3,1)=lat(3,1)-round(lat(3,1)/lat(1,1))*lat(1,1);
coor=pos*(lat^-1);
coor=coor-floor(coor);

%write coo file
fprintf(fp, '0.0 %6.7f xlo xhi\n', lat(1,1)+1e-5);
fprintf(fp, '0.0 %6.7f ylo yhi\n', lat(2,2)+1e-5);
fprintf(fp, '0.0 %6.7f zlo zhi\n', lat(3,3)+1e-5);
fprintf(fp, '%6.7f %6.7f %6.7f xy xz yz\n',lat(2,1),lat(3,1),lat(3,2));
fprintf(fp, '\n');

fprintf(fp, 'Atoms\n');
fprintf(fp, '\n');
coordLoop = 1;
pos=coor*lat;
for i = 1 : length(numIons)
  for j = 1 : numIons(i)
    fprintf(fp, '%d  %d %8.4f %8.4f %8.4f\n', coordLoop, i,pos(coordLoop,:));
    coordLoop = coordLoop + 1;
  end
end
fprintf(fp, '\n');


fclose(fp);
[nothing, nothing] = unix('rm out.xyz lammps.out');
