function vol = calcDefaultVolume(numIons, atomType, targetPressure, isMol, tempFolder, codeFolder)


if ~exist(tempFolder,'dir')
    [nothing, nothing] = unix(['mkdir -p '  tempFolder])	;
end
fp = fopen([tempFolder '/vol.input'],'w');
fprintf(fp, '%5d\n', length(numIons));
fprintf(fp, '%5d ', numIons);
fprintf(fp, '\n ');
fprintf(fp, '%5d ', atomType);
fprintf(fp, '\n ');
fprintf(fp, '%8.3f\n', targetPressure);
fprintf(fp, '%2d\n', isMol);
fclose(fp);

cd(tempFolder);

[nothing, result] = unix([codeFolder '/calcDefaultVolume < ' tempFolder '/vol.input' ]);

if ~isempty( result )
    vol = str2num(result);
else
    disp(' Error in calling calcDefaultVolume, exit USPEX ...')
    exit
end

cd ..
