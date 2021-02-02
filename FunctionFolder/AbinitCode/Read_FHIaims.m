function target = Read_FHIaims(flag, ID)
if     flag == -1
[nothing, results] = unix('grep "Have a nice day" FHI_output');
if isempty(results)
disp('FHI-aims is not completely Done');
[nothing, nothing] = unix(['cp FHI_output ERROR-OUTPUT-' ID]);
target = 0;
else
target = 1;
end
elseif flag == 0
[nothing, results] = unix('./getStuff FHI_output "Total energy corrected" 7');
if isempty(results)
disp('FHI-aims is not done correctly!!!'); 
USPEXmessage(1351,'',0);
[nothing, nothing] = unix(['cp FHI_output ERROR-OUTPUT-' ID]);
target = 0;
elseif ~exist('FHI_output')
disp('FHI Output does not exist!');
target = 0 ;
else
target = 1;
end
elseif flag ==  1 
[nothing, results] = unix('./getStuff FHI_output "Total energy corrected" 7 | tail -1'); 
target = str2num(results);
end
