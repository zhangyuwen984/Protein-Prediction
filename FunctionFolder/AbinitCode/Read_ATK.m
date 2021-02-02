function target = Read_ATK(flag)
%-1: if ATK is complete
%0: if 1st SCF is done
% 1: energy/enthalpy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag == -1  
	  [nothing, results] = unix('grep Timing ATK.out');
	  if isempty(results)
	    disp('ATK calculation is not completely Done');
	    target = 0;
	  else
	    target = 1;
	  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
elseif flag == 0
	[nothing, results] = unix('./getStuff ATK.out "Total energy " 5');
	if isempty(results)
		disp('PROBLEM: can not read the energy from ATK.out for structure number ');
		target = 0;
	else
		target = 1;
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  1 
	[nothing, Energy_Str] = unix('./getStuff ATK.out "Total energy " 5');
	Energy = str2num(Energy_Str);
	target = Energy(end);
end
