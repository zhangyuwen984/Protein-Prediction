function [coor, lat] = Read_Structure(code, Const_Lat)
%dimension: -3 -2 0 1 2 3
%molecule: 0/1
if code == 1
   [coor, lat] =    Read_VASP_Structure(); 
elseif code ==2
   [coor, lat] =  Read_SIESTA_Structure();
elseif code ==3
   [coor, lat] =    Read_GULP_Structure();
elseif code ==4
   [coor, lat] =  Read_LAMMPS_Structure(); 
elseif code ==5
   [coor, lat] =Read_NeurNetw_Structure();
elseif code ==6
   [coor, lat] = Read_DMACRYS_Structure();
elseif code ==7
   [coor, lat] =    Read_CP2K_Structure();
elseif code ==8
   [coor, lat] =      Read_QE_Structure();
elseif code ==9
   [coor, lat] = Read_FHIaims_Structure();
elseif code ==10
   [coor, lat] =     Read_ATK_Structure();
elseif code ==11
   [coor, lat] =  Read_CASTEP_Structure();
elseif code == 12
   [coor, lat] =  Read_Tinker_Structure_MOL();
elseif code == 13
   [coor, lat] =  Read_MOPAC_Structure();
elseif code == 14
   [coor, lat] =    Read_VASP_Structure();
end

if (Const_Lat==0) 
   if code~=2 && code~=12 % SIESTA and Tinker
       coor = coor*lat;
      [coor,lat] = optLattice(coor, lat);
       coor = coor/lat;
   end
end
