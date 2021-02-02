function stillRunningInfo()


tmp = dir('still_running');
if exist('still_running')
    if (now - tmp.datenum > 1/12)    % 1.00 = 24 hours
        [nothing, nothing] = unix('rm still_running');
    else
        disp('still_runing is present, matlab has to exit, please see NOT_YET file for details');
        [nothing, nothing] = unix('echo "This is a warning message"   >> NOT_YET');
        [nothing, nothing] = unix('echo "It indicates that you attempted to call matlab"      >> NOT_YET');
        [nothing, nothing] = unix('echo "when file still_running is still present"            >> NOT_YET');
        [nothing, nothing] = unix('echo "It would be problematic when the numParallel is on"  >> NOT_YET');
        [nothing, nothing] = unix('echo "Possible reasons: "                                  >> NOT_YET');
        [nothing, nothing] = unix('echo "1, The time interval to call MATLAB is too short"    >> NOT_YET');
        [nothing, nothing] = unix('echo "2, MATLAB exits with error"                          >> NOT_YET');
        quit
    end
end

[nothing, nothing] = unix('echo "This file is present for two possible reasons"      > still_running');
[nothing, nothing] = unix('echo "1, MATLAB is still running"                        >> still_running');
[nothing, nothing] = unix('echo "2, MATLAB exits with error"                        >> still_running');
[nothing, nothing] = unix('echo "If it stays for long time when numParallel is on"  >> still_running');
[nothing, nothing] = unix('echo "Matlab is either in dead loop or exits with error.">> still_running');
[nothing, nothing] = unix('echo "please do check it                              "  >> still_running');
