function USPEXmessage(msgID, msgString, option)

%  This function is used to display the warnings on screen or write into the files
%
%      msgID : 1) The warning message ID
%              2) if msgID = 0, display the contest of msgString
%              3) if msgID > 0, display the default message
%
%  msgString : the warning message defined by user
%
%     option : display option for output of warning message
%              1) if option = 0, display message on screen and write it in file
%              2) if option < 0, display message on screen
%              3) if option > 0, write message in file
%  added by Guangrui Qian

global ORG_STRUC


clockTIME=clock;
timeString = ['===== USPEX WARNING @ ',num2str(clockTIME(4)),':',num2str(clockTIME(5)),':',num2str(clockTIME(3)), ' ', date];

if msgID > 0	
    displayMsg = [ timeString, ' ===== \n', messageInfo(msgID)];
else
	displayMsg = [ timeString, ' ===== \n'];
end

if length(msgString)>0
    displayMsg = [displayMsg, msgString, '\n\n'];
else
    displayMsg = [displayMsg, '\n\n'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display the warning message on screen
if option <= 0
    fprintf(1, displayMsg);
end
%% write the warning message in file
if option >=0
    filePath=[ORG_STRUC.homePath, '/Warnings'];
    fp = fopen(filePath, 'a+');
    fprintf(fp, displayMsg);
    fclose(fp);
end

