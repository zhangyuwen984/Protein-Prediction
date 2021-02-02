function target = Read_CP2K(flag)

%-1: if CP2K is complete
% 0: if 1st SCF is done & SCF is converged
% 1: energy/enthalpy


Ha = 27.211385; % eV


if     flag == -1
    if exist('USPEX-pos-1.xyz','file') && exist('USPEX-1.cell','file')
        target = 1;
    else
        target = 0;
        messageInfo(1301);
    end
elseif flag == 0
    try
        [nothing, cellString] = unix('tail -n 1 USPEX-1.cell');
        lat = str2num(cellString);
        if length(lat)==12
            target = 1;
        else
            target = 0;
            messageInfo(1302);
        end
    catch
        target = 0;
        messageInfo(1302);
    end
elseif flag ==  1 %
    [nothing, Energy_Str] = unix(['./getStuff USPEX-pos-1.xyz "E =" 10 ']);
    Energy = str2num(Energy_Str);
    if isempty(Energy)
        [nothing, Energy_Str] = unix(['./getStuff USPEX-pos-1.xyz "E =" 7 ']);
    end
    
    if length (Energy_Str) > 21
        Energy_Str = Energy_Str(end-20:end);
    end
    Energy = str2num(Energy_Str);
    target = Energy(end)*Ha; % eV;
end
