function [coor, lat] = Read_MOPAC_Structure()
global ORG_STRUC
coor = callAWK('MOPAC_coor.awk', 'input.arc');
if ORG_STRUC.dimension == 0
   lat = [30 0 0; 0 30 0; 0 0 30];
   Xcm=(sum(coor(:,1)))/length(coor(:,1));
   Ycm=(sum(coor(:,2)))/length(coor(:,2));
   Zcm=(sum(coor(:,3)))/length(coor(:,3));
   Dif= [15-Xcm 15-Ycm 15-Zcm];
   for i = 1 : 3
        coor(:,i) = coor(:,i) + Dif(i);
   end
elseif ORG_STRUC.dimension == 3
   lat = callAWK('MOPAC_lat.awk', 'input.arc');
end
coor = coor/lat;
