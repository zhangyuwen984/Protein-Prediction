function Write_AbinitCode(code, Ind)
if code == 1
    Write_VASP(Ind)
elseif code ==2
    Write_SIESTA(Ind)
elseif code ==3
    Write_GULP(Ind)
elseif code ==4
    Write_LAMMPS(Ind)
elseif code ==5
    Write_NeurNetw(Ind)
elseif code ==6
    Write_DMACRYS(Ind)
elseif code ==7
    Write_CP2K(Ind)
elseif code ==8
    Write_QE(Ind)
elseif code ==9
    Write_FHIaims(Ind)
elseif code ==10
    Write_ATK(Ind)
elseif code ==11
    Write_CASTEP(Ind)
elseif code == 12
    Write_Tinker(Ind)
elseif code == 13
    Write_MOPAC(Ind)
elseif code == 14
    Write_BoltzTraP(Ind)
end


