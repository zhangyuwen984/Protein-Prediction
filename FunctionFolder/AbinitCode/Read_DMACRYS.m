function target = Read_DMACRYS(flag)
%-1: if VASP is complete
%-0: if 1st SCF is done
% 1: energy/enthalpy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if     flag == -1
       [nothing, results] = unix('grep "Total run time" output');
       if isempty(results)
          disp('DMACRYS is not completely Done');
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag == 0
       [nothing, results] = unix(['./getStuff fort.12 Energy 4']);
       if isempty(results) | isempty(str2num(results))
          disp('DMACRYS 1st SCF is not done');
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  1 %
       [nothing, Energy_Str] = unix(['./getStuff fort.12 Energy 4 |tail -1']);
       target = str2num(Energy_Str);
end
