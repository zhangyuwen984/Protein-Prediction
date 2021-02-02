function writeMagmoment(resFolder)

global ORG_STRUC
global USPEX_STRUC

if ORG_STRUC.spin > 0

   fpath = [ resFolder '/Individuals_magmoment'];
   %fpath2= [ resFolder '/magmoment_complete.dat'];
   fp = fopen(fpath, 'w');
   %fp2= fopen(fpath2,'w');
   if (ORG_STRUC.maxAt==ORG_STRUC.minAt) | (ORG_STRUC.varcomp == 0)
       elementsStr = createElementsStr(ORG_STRUC.numIons,ORG_STRUC.atomType);
   else
       elementsStr =[];
   end
   fprintf(fp, ['  ID  MagType    -Atomic Magmoments-\n']);
   fprintf(fp, ['               ' elementsStr '\n']);
   for i = 1:length( USPEX_STRUC.POPULATION )
        magTypeIni = magTypeString( USPEX_STRUC.POPULATION(i).magmom_ions(1,1) );
        magTypeFin = magTypeString( USPEX_STRUC.POPULATION(i).magmom_ions(end,1) );
        fprintf(fp, '%4d %6s : ', i, magTypeFin );
        fprintf(fp, '%8.3f ', USPEX_STRUC.POPULATION(i).magmom_ions(end, 2:end));
        fprintf(fp, '\n');
   %=============================================================================
   %     for 
   %     fprintf(fp, '%4d %6s->%6s: ', i, magTypeIni, magTypeFin );
   %     fprintf(fp, '%8.3f ', USPEX_STRUC.POPULATION(i).magmom_ions(end, 2:end));
   %     fprintf(fp, '\n');
   end
   fclose(fp);
   %fclose(fp2);
end


%-----------------------------------------------------------
function elementsStr = createElementsStr(numIons,atomType)

elementsStr=[];
for i = 1:length(numIons)
    for j = 1:numIons(i)
        elementsStr=[elementsStr '    ' megaDoof(atomType(i)) '    '];
    end
end


