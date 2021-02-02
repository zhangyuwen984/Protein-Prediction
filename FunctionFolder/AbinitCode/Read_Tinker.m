function target = Read_Tinker(flag)
%-0: if 1st SCF is done
% 1: energy/enthalpy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('output.frac')
   if flag == 0
       [nothing, results] = unix(['grep "Function Value"  output']);
       if isempty(results) 
          disp('TINKER  1st SCF is not done');
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   elseif flag ==  1 %
       [nothing, Energy_Str] = unix(['grep "Function Value"  output |tail -1']);
       target = str2num(Energy_Str(end-10:end));
       if target < -1e5
          target=10000;
       end
   end
else
   if flag == 0
          [nothing, results] = unix(['./getStuff output "The lowest energy:" 4']);
          if isempty(results) || isempty(str2num(results))
             disp('Tinker 1st SCF is not done');
             target = 0;
          else
             target = 1;
          end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   elseif flag ==  1 %
          [nothing, Energy_Str] = unix(['./getStuff output "The lowest energy:" 4']);
          target = str2num(Energy_Str);
          if target < -1e5
              target=10000;
          end

   end
end
