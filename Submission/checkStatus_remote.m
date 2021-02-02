function doneOr = checkStatus_remote(jobID, USPEX, Folder)

global POP_STRUC
global ORG_STRUC

%--------------------------------------------------------------------
%This routine is to check if the submitted job is done or not
%One needs to do a little edit based on your own case.
%--------------------------------------------------------------------

%Step1: Specify the PATH to put your calculation folder
Home = ['/home/AI/shipilov.ab/statistic_collection/1cei']; %'pwd' of your home directory of your remote machine
Address = 'shipilov.ab@calc.cod.phystech.edu'; %your target server: username@address
Path = [Home '/' USPEX '/CalcFold' num2str(Folder)];

%Step2: Check JobID, the exact command to check job by jobID
[nothing,nothing]=unix(['ssh ' Address ' "squeue" > ' ORG_STRUC.homePath '/jobinfo.dat']);

%    [nothing, nothing] = unix(['squeue > '  ORG_STRUC.homePath '/jobinfo.dat']);
    [nothing, statusStr] = unix(['cat ' ORG_STRUC.homePath '/jobinfo.dat | grep -v CG | grep -v JOBID | grep ' num2str(jobID) ' | wc -l '] );
    [nothing, nothing] = unix(['cat ' ORG_STRUC.homePath '/jobinfo.dat | grep -v CG | grep -v JOBID | grep ' num2str(jobID) ' '] );
    % disp(statusStr)
    isOK=str2num( statusStr );

    if  isOK==0
        doneOr = 1;
        [nothing, nothing] = unix(['scp -r ' Address ':' Path '/* .']);
        [nothing, nothing] = unix(['ssh ' Address ' "rm -r" ' Path ]);
    else
        doneOr = 0;
    end

