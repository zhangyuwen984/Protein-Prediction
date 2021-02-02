function target = Read_SIESTA(flag)
%-1: if SIESTA is complete
% 0: if 1st SCF is done
% 1: energy/enthalpy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if     flag == -1
       [nothing, results] = unix('grep "End of run" output');
       if isempty(results)
          disp('SIESTA output is not completely Done');
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag == 0
      [nothing, results] = unix(['./getStuff output enthalpy 4']);
       if isempty(results)
          disp('SIESTA 1st SCF is not done');
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
elseif flag ==  1 %
      [nothing, Energy_Str] = unix(['./getStuff output enthalpy 4']);
       Energy = str2num(Energy_Str);
       target = Energy(end);
end
