%This function return a number as the value of property  
% Example: Seebeck coefficient 
% Dong Dong 
%

% function [property] = calcProperty()
function [property] = calcThermoelectricProperties()
    Temp0 = 300.0;
    property = -1.0;

    % detemine the current folder name
    [upperpath,directory] = fileparts(pwd());

    % determine Seebeck Coefficient from *.trace
    trace_files = dir([directory '-*.trace']);
    fp=fopen(trace_files(1).name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max seebeck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SeebeckF(1) = 0.0;
% tmp = fgetl(fp); 
% i = 1;
% while ~feof(fp)
%     tmp=fgetl(fp);
%     tmp=str2num(tmp);
%     if ( tmp(2) > (Temp0-5.0) ) & (tmp(2) < (Temp0+5.0) )
%         Temp(i) = tmp(2);
%         SeebeckF(i) = tmp(5);
%         i=i+1;
%     end
% end
% fclose(fp);

 % read the maximum of seebeck in muV/K
% property = 1000000 * max(SeebeckF); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seebeck at the fermi level (mu=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power factor at the fermi level (mu=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SeebeckF(1) = 0.0;
SeebeckF(2) = 0.0;
SigmaF(1) = 0.0;
SigmaF(2) = 0.0;
mu(1) = -10.0;
mu(2) = -10.0; 

tmp = fgetl(fp); 
while ~feof(fp)
    tmp=fgetl(fp);
    tmp=str2num(tmp);
    if ( tmp(2) > (Temp0-5.0) ) & (tmp(2) < (Temp0+5.0) )
        mu(1) = tmp(1); 
        SeebeckF(1) = tmp(5);
        SigmaF(1) = tmp(6);
    end
    if ( mu(2) < 0.0 ) & ( mu(1) < 0.0 )
        mu(2)= mu(1);
        SeebeckF(2)=SeebeckF(1);
        SigmaF(2) = SigmaF(1);
    else 
        break;
    end
end

seebeck = (SeebeckF(2)*mu(1)-SeebeckF(1)*mu(2))/(mu(1)-mu(2)); 
sigma = (SigmaF(2)*mu(1)-SigmaF(1)*mu(2))/(mu(1)-mu(2));

% seebeck
% property = 1000000 * seebeck; 

% power factor P/tau*10e10
property = 10e-10*seebeck*seebeck*sigma;

fclose(fp);

