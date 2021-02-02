function WriteIndividual(resFolder, energyUnits)

%To update all the necessary items at the end of each generation
%Individuals would be updated
%Last updated by MR (2014/04/10)
global ORG_STRUC
global USPEX_STRUC

if nargin == 2 && strcmp(energyUnits, 'kcal/mol')
    energyUnits = '(kcal/mol)';
else
    energyUnits     = '   (eV)   ';
end


fpath = [ resFolder '/Individuals'];
fp = fopen(fpath, 'w');
if ORG_STRUC.spin == 1
    fprintf(fp,  'Gen   ID    Origin   Composition    Enthalpy       RMSD     Fitness   Q_entr A_order S_order Magmom-Type\n');
else
    fprintf(fp,  'Gen   ID    Origin   Compositios    Enthalpy       RMSD     Fitness    Q_entr A_order S_order\n');
end
fprintf(fp, ['                                   ' energyUnits '      nm\n']);

for i=1:length(USPEX_STRUC.POPULATION)
      gen    = USPEX_STRUC.POPULATION(i).gen;
    enth     = USPEX_STRUC.POPULATION(i).Enthalpies(end);
    num      = USPEX_STRUC.POPULATION(i).numIons;
    fit      = USPEX_STRUC.POPULATION(i).Fitness;
    howcome  = USPEX_STRUC.POPULATION(i).howCome;
    entropy  = USPEX_STRUC.POPULATION(i).struc_entr;
          s  = USPEX_STRUC.POPULATION(i).S_order;
      order  = USPEX_STRUC.POPULATION(i).order;
    rmsd     = USPEX_STRUC.POPULATION(i).RMSD;
    if ORG_STRUC.spin == 1
        magmom   = sum(USPEX_STRUC.POPULATION(i).magmom_ions(end,2:end));
        magType  = USPEX_STRUC.POPULATION(i).magmom_ions(end,1);
    else
        magmom   =[];
        magType  =[];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process 'N/A' case:
    comp_format = '[%11s]';
    if strcmp(num, 'N/A')
        composition = num;
        comp_format = '%-13s';
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

    if isempty(entropy)
       entropy = 0;
    end
    
    if strcmp(order, 'N/A')
        a_o = order;
    else
        if sum(num)>0
           a_o = sum(order)/sum(num);
        else
           a_o = 0;
        end
    end

    fprintf(fp,['%3d %4d %-11s '  comp_format ' %9.3f ' ' %9.3f  '  ' %10.3f '  ' ' entropy_format ' ' ao_format ' ' s_format '\n'], ...
                 gen, i, howcome, composition,  enth,      rmsd,       fit,           entropy,           a_o,          s);
end
fclose(fp);

