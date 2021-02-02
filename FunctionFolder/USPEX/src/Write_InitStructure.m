function Write_InitStructure(Count, INIT_numIons, INIT_LAT, INIT_COORD, resFolder)

%Output the structure info before relaxation
%Lastly updated by Qiang Zhu (2014/02/18)

fpath3 = [ resFolder '/gatheredPOSCARS_unrelaxed' ];
fp3=fopen(fpath3, 'a+');

fprintf(fp3, 'Structure %4d\n', Count);
fprintf(fp3, '1.000000\n');

for latticeLoop = 1 : 3
    fprintf(fp3, '%12.6f %12.6f %12.6f\n', INIT_LAT(latticeLoop,:));
end
for i=1:length(INIT_numIons)
    fprintf(fp3, '%4d ', INIT_numIons(i));
end
fprintf(fp3, '\n');
fprintf(fp3, 'Direct\n');

for coordLoop = 1 : sum(INIT_numIons)
 fprintf(fp3, '%12.6f %12.6f %12.6f\n', INIT_COORD(coordLoop,:));
end
fclose(fp3);
