function target = Read_QE(flag, ID)
%This rountine is to read output from QE
%File: output (4.6/5.2)
%-1: QE is complete (FOR NEB)
% 0: 1st SCF is done (For USPEX, etc)
% 1: Energy 
% 2: pressue tensor 
% 7: atomic forces
%Last updated by Qiang Zhu (2013/10/03)

Ry = 13.6056923; % eV
bohr= 0.529177249;%Angstron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if     flag == -1  % to change!!!
       [nothing, results] = unix('grep "JOB DONE" output');
       if isempty(results)
          disp('Quantum Espresso is not completely Done');
          [nothing, nothing] = unix(['cp output ERROR-OUTPUT-' ID]);
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag == 0
          [nothing, results1] = unix('./getStuff output "enthalpy new" 5'); % enthalpy in Ry
          [nothing, results2] = unix('./getStuff output "!    total energy" 5'); % energy in Ry, 
       if isempty(results1) & isempty(results2)
          disp('Quantum Espresso 1st SCF is not done');
          [nothing, nothing] = unix(['cp output ERROR-OUTPUT-' ID]);
          target = 0;
       else
          target = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
elseif flag ==  1 %
        [nothing, results1] = unix('./getStuff output "enthalpy new" 5 |tail -1'); % enthalpy in Ry
       if ~isempty(results1) 
           target = str2num(results1);
       else
          [nothing, results2] = unix('./getStuff output "!    total energy" 5 |tail -1'); % energy in Ry
           target = str2num(results2);
       end
       target = target*Ry;
elseif flag ==  2 %
       target = callAWK('QE_pres.awk','output');
elseif flag ==  7%
    try
        [nothing, numStr] = unix(['grep "atoms/cell" output |cut -d"=" -f2']);
        numIons= sum( str2num(numStr) );
        force_orig = callAWK('QE_force.awk', 'output', ['num=', num2str(numIons)]);
        if isempty(force_orig)
            %USPEXmessage(1109,'',0);
            target = [];
        else
            target = force_orig*Ry/bohr;
        end
    catch
        target = [];
    end
end
