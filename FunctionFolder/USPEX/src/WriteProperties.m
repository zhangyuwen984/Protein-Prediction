function WriteProperties(type, resFolder)
%Properties output
%type: 3 hardness; 6 dielectric; 7 bandgap; 8 d+gap; 9: magmoment; 1101-1109: elastic properties
%Last updated by Qiang Zhu (2014/02/18)
global USPEX_STRUC

fpath = [ resFolder '/Properties'];
fp = fopen(fpath, 'w');
if type == 3 %hardness
   fprintf(fp, 'hardness\n');
   fprintf(fp, ' (GPa)  \n');
   for i=1:length(USPEX_STRUC.POPULATION)
       fprintf(fp, '%8.3f\n', USPEX_STRUC.POPULATION(i).hardness);
   end
elseif type == 6 %dielectric
   fprintf(fp, '                    dielectric tensor                      trace    \n');
   fprintf(fp, '     11       22       33       12       13       23                \n');
   for i=1:length(USPEX_STRUC.POPULATION)
       tensor = USPEX_STRUC.POPULATION(i).dielectric_tensor;
	   trace = sum(tensor(1:3))/3;
       fprintf(fp, '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', tensor(:), trace);
   end
elseif type==7 %bandgap
   fprintf(fp, 'band_gap\n');
   fprintf(fp, ' (eV)\n');
   for i=1:length(USPEX_STRUC.POPULATION)
       fprintf(fp, '%6.3f\n', USPEX_STRUC.POPULATION(i).gap);
   end
elseif type==8 %Tr(dielectric tensor)/3 multiplied by gap^2 value
   fprintf(fp, '    E_g    lamda   lamda^3*(E_g/4)^2\n');
   Egc = 4;
   for i=1:length(USPEX_STRUC.POPULATION)
       gap   = USPEX_STRUC.POPULATION(i).gap;
       lamda = sum(USPEX_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
       if gap >= Egc
          fit = (lamda^3)*(gap/Egc)^2;
       else
          fit = (lamda^3)*(gap/Egc)^6;
       end
       fprintf(fp, '%8.3f %8.3f %8.3f\n', gap, lamda, fit);
   end
elseif type==9
   fprintf(fp, 'MagMoment\n');
   fprintf(fp, ' (eV)\n');
   for i=1:length(USPEX_STRUC.POPULATION)
       fprintf(fp, '%6.3f\n', USPEX_STRUC.POPULATION(i).gap);
   end
elseif (type>1100) & (type<1112)  % elastic related properties
   fpath0 = [ resFolder '/AuxiliaryFiles/ElasticConstantMatrix'];
   fp0 = fopen(fpath0, 'w');
   
   fprintf(fp, 'ID   Modulus:Bulk, Shear, Youngs  Ratio:Poissons,Pughs Vicker-Hard Toughness Debye-T Velocity:Mean,S-wave,P-wave Stable\n');
   fprintf(fp, '            (GPa)  (GPa)   (GPa)                           (GPa)  (MPa*m^1/2)  (K)       (km/s)  (km/s)  (km/s)    \n');
   
   fprintf(fp0,'Elastic Constant Matrix (GPa)\n');
   for i = 1:length(USPEX_STRUC.POPULATION)
	   fprintf(fp, ' %4d    %8.1f%7.1f%8.1f        %6.3f  %6.3f  %8.2f  %8.2f  %7.2f   %8.2f %7.2f %7.2f ', i, USPEX_STRUC.POPULATION(i).elasticProperties(1:end-1));
	   if USPEX_STRUC.POPULATION(i).elasticProperties(end) == 1
		   fprintf(fp, '   .Y. \n');
	   else
		   fprintf(fp, '   .N. \n');
	   end
	   
	   fprintf(fp0,'EA%d\n', i);
	   if isempty(USPEX_STRUC.POPULATION(i).elasticMatrix)
		   fprintf(fp0,'   NaN: no availiable data\n');
	   else
	   	   %USPEX_STRUC.POPULATION(i).elasticMatrix
	   	   for j = 1:6
	           fprintf(fp0,'   %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n', USPEX_STRUC.POPULATION(i).elasticMatrix(j,:));
		   end
	   end
   end	
   fclose(fp0);   
   %%%%% print the elastic tensor
else
   fprintf(fp, 'See Individuals\n');
end
fclose(fp);



%%---------
if exist('writeMagmoment')
    writeMagmoment(resFolder);
end
