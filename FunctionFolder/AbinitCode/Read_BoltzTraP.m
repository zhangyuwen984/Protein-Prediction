%This function check if the property calculation succeeds or not 
% Success: 1; No: 0;
% Dong Dong 

% function [target] = checkProperty()
function [target] = Read_BoltzTraP(flag, ID)

    % detemine the current folder name
    currentdir = pwd;
    [upperpath,directory] = fileparts(currentdir);
  
    % Added by Fei Qi, 2015/07/27
    if flag == 0
        trace_files = dir([directory '-*.trace']);
        if length(trace_files) >= 1
            target = 1;
        else
            target = 0;
        end
    elseif flag == 2
        target = calcThermoelectricProperties();
    end

  % Commented out Dong Dong's implementation. %% Fei Qi, 2015/07/27
  % if exist([dir '.trace'],'file')
  %    target = 1;
  % else
  %    target = 0;
  % end

