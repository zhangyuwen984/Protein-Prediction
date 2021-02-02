function plotE(resFolder)
global USPEX_STRUC
global ORG_STRUC

cd (resFolder)
N = length(USPEX_STRUC.GENERATION);
E = zeros(N,1);
for i = 1:N
   E(i) = USPEX_STRUC.GENERATION(i).Fitness;
end

numImage=length(E);
MinE=min(E);
MaxE=max(E);
dE=(MaxE-MinE)/(numImage);
if dE==0
   dE = 0.001;
end
Erescale=round( (E-MinE)/dE );

Nmatrix=numImage;
for i=1:numImage
   LINE = Nmatrix+1-Erescale(i);
   FIG(LINE,(i-1)*3+3)=1;
end
FIG(:,2)=2;
%--------------------------------------------------------

fpath=['FIGURE'];
fp=fopen(fpath,'a+');
for i=1:Nmatrix+1
    fprintf(fp,'%1d',FIG(i,:));
    fprintf(fp,'\n');
end
fclose(fp);
[nothing, nothing] = unix( ['sed -i "s/0/ /g" FIGURE'] );
[nothing, nothing] = unix( ['sed -i "s/1/+/g" FIGURE'] );
[nothing, nothing] = unix( ['sed -i "s/2/|/g" FIGURE'] );

fp=fopen(fpath,'a');
for i=1:numImage
 fprintf(fp,'%3s','--o');
end
fprintf(fp,'\n');
for i=1:numImage
 fprintf(fp,'%3d',i);
end
fprintf(fp,'\n\n');

if ORG_STRUC.dimension == -4 && ORG_STRUC.varcomp == 0 && ORG_STRUC.molecule == 0
    units = 'kcal/mol';
else
    units = 'eV';
end

fprintf(fp,'The lowest target fitness is %6.4f %s\n',E(numImage), units);

fclose(fp);
[nothing, nothing] = unix( ['cat FIGURE >> OUTPUT.txt'] );
[nothing, nothing] = unix( ['rm FIGURE'] );

cd ..
%---------------------------------------------------------
%--     FUNCTION END
