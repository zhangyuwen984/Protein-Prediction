function WriteOUTPUT(Ind_No, resFolder)

%To update all the necessary items after each structure is done
%1, OUTPUT.txt
%2, Individuals
%3, Origin
%Lastly updated by Qiang Zhu (2014/02/18)
global ORG_STRUC
global USPEX_STRUC

  gen    = USPEX_STRUC.POPULATION(Ind_No).gen;
symg     = USPEX_STRUC.POPULATION(Ind_No).symg;
enth     = USPEX_STRUC.POPULATION(Ind_No).Enthalpies(end);
num      = USPEX_STRUC.POPULATION(Ind_No).numIons;
%fit      = USPEX_STRUC.POPULATION(Ind_No).Fitness;
howcome  = USPEX_STRUC.POPULATION(Ind_No).howCome;
KPOINTS  = USPEX_STRUC.POPULATION(Ind_No).K_POINTS(end,:);
volume   = USPEX_STRUC.POPULATION(Ind_No).Vol;
density  = USPEX_STRUC.POPULATION(Ind_No).density;
 entropy = USPEX_STRUC.POPULATION(Ind_No).struc_entr;
       s = USPEX_STRUC.POPULATION(Ind_No).S_order;
   order = USPEX_STRUC.POPULATION(Ind_No).order;
if ORG_STRUC.spin == 1
    magmom   = sum(USPEX_STRUC.POPULATION(Ind_No).magmom_ions(end,2:end));
    magType  = USPEX_STRUC.POPULATION(Ind_No).magmom_ions(end,1);
else
    magmom   = [];
    magType  = [];
end
sum_num = sum(num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process 'N/A' case:
comp_format = '[%11s]';
if strcmp(num, 'N/A')
    composition = num;
    sum_num = 1;
    comp_format = '%-13s';
end

volume_format = '%9.3f';
if strcmp(volume, 'N/A')
    volume_format = '  %-7s';
end

density_format = '%7.3f';
if strcmp(density, 'N/A')
    density_format = '%-6s';
end

kpoints_format = '[%2d %2d %2d]';
if strcmp(KPOINTS, 'N/A')
    kpoints_format = '  %-8s';
end

symmetry_format = '%3d';
if strcmp(symg, 'N/A')
    symmetry_format = '%-5s';
end

entropy_format = '%6.3f';
if strcmp(entropy, 'N/A')
    entropy_format = '%-6s';
end

ao_format = '%6.3f';
if strcmp(order, 'N/A')
    ao_format = '%-7s';
end

s_format = '%6.3f';
if strcmp(s, 'N/A')
    s_format = '%-6s';
end


    mag_format = '%8.3f';
    if strcmp(magmom, 'N/A')
        mag_format = '%-6s';
    end

    magType_format = '%6s';
    if strcmp(magmom, 'N/A')
        magType_format = '%-6s';
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(num, 'N/A')
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
       composition=[composition,blanks(shift(length(num)))];
    end
end

if ~isempty(USPEX_STRUC.POPULATION(Ind_No).Parents)
   par_ID = USPEX_STRUC.POPULATION(Ind_No).Parents.parent;
   par_fit= USPEX_STRUC.POPULATION(Ind_No).Parents.enthalpy;
else
   par_ID = '0';
   par_fit = enth/sum_num;
end

if isempty(entropy)
   entropy = 0;
end

if strcmp(order, 'N/A')
    a_o = order;
else
    if sum_num>0
       a_o = sum(order)/sum_num;
    else
       a_o = 0;
    end
end

fpath  = [resFolder '/OUTPUT.txt'];
fpath1 = [resFolder '/Individuals'];
fpath2 = [resFolder '/origin'];
fp  = fopen(fpath,  'a+');
fp1 = fopen(fpath1, 'a+');
fp2 = fopen(fpath2, 'a+');

fprintf(fp,['%4d %-11s ' comp_format ' %10.3f  ' volume_format '   ' kpoints_format '  ' symmetry_format '\n'], Ind_No, howcome, composition, enth, volume, KPOINTS(:), symg);

if ORG_STRUC.spin == 1
fprintf(fp1, ['%3d %4d %-11s '      comp_format ' %9.3f ' volume_format ' ' density_format '     N/A    ' kpoints_format ' ' symmetry_format ' ' entropy_format ' ' ao_format ' ' s_format ' ' mag_format ' ' magType_format '\n'], ...
              gen, Ind_No, howcome, composition, enth,    volume,           density,                      KPOINTS(:),        symg,               entropy,           a_o,          s, magmom, magTypeString(magType(1)));
else
fprintf(fp1, ['%3d %4d %-11s '      comp_format ' %9.3f ' volume_format ' ' density_format '     N/A    ' kpoints_format ' ' symmetry_format ' ' entropy_format ' ' ao_format ' ' s_format '\n'], ...
              gen, Ind_No, howcome, composition, enth,    volume,           density,                      KPOINTS(:),        symg,               entropy,           a_o,          s);

end

fprintf(fp2,'%4d %-11s %8.3f  %8.3f  [%10s]\n', Ind_No, howcome, enth/sum_num, par_fit, par_ID);
fclose(fp);
fclose(fp1);
fclose(fp2);

[nothing, nothing] = unix(['echo ' num2str(USPEX_STRUC.POPULATION(Ind_No).Enthalpies,'%10.3f') ' >> ' resFolder '/enthalpies_complete.dat']);
