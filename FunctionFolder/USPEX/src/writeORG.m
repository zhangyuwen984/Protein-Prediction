function writeORG()

global ORG_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

if ORG_STRUC.molecule == 1
    fprintf(fp,[alignLine('Molecular Crystals') '\n']);
    text = {'Zhu Q., Oganov A.R., Glass C.W., Stokes H. (2012)', ...
        'Constrained evolutionary algorithm for structure prediction of', ...
        'molecular crystals: methodology and applications.', ...
        'Acta Cryst. B, 68, 215-226' ...
        };
    formatted_rows = createHeader_wrap(text, 'left');
    for i=1:size(formatted_rows,2)
        fprintf(fp,[formatted_rows{i} '\n']);
    end
end

if ORG_STRUC.dimension == 2
    fprintf(fp,[alignLine('Surfaces') '\n']);
    text = {'Zhu Q., Li L., Oganov A.R., Allen P.B. (2013)', ...
        'Evolutionary method for predicting surface reconstructions', ...
        'with variable stoichiometry.', ...
        'Phys. Rev. B, 87, 195317' ...
        };
    formatted_rows = createHeader_wrap(text, 'left');
    for i=1:size(formatted_rows,2)
        fprintf(fp,[formatted_rows{i} '\n']);
    end
end

if ORG_STRUC.dimension == 1
    fprintf(fp,[alignLine('Polymers') '\n']);
    text = {'Zhu, Q., Sharma V., Oganov A.R., Ramprasad R. (2014)', ...
        'Predicting polymeric crystal structures by evolutionary algorithms.', ...
        'J. Chem. Phys. 141, 154102' ...
        };
    formatted_rows = createHeader_wrap(text, 'left');
    for i=1:size(formatted_rows,2)
        fprintf(fp,[formatted_rows{i} '\n']);
    end
end

if ORG_STRUC.dimension == -2
    fprintf(fp,[alignLine('2D-crystals') '\n']);
    text = {'Zhou X.-F., Dong X., Oganov A.R., Zhu Q.,', ...
        'Tian Y., Wang H.-T. (2014)', ...
        'Semimetallic Two-Dimensional Boron Allotrope', ...
        'with Massless Dirac Fermions.', ...
        'Phys. Rev. Lett. 112, 085502' ...
        };
    formatted_rows = createHeader_wrap(text, 'left');
    for i=1:size(formatted_rows,2)
        fprintf(fp,[formatted_rows{i} '\n']);
    end
end


if ORG_STRUC.varcomp == 1
    fprintf(fp,[alignLine('Variable Composition') '\n']);
    text = {'Lyakhov A.O., Oganov A.R., Valle M. (2010)', ...
        'Crystal structure prediction using evolutionary approach.', ...
        'In: Modern methods of crystal structure prediction (ed: A.R. Oganov)', ...
        'Berlin: Wiley-VCH', ...
        '', ...
        'Oganov A.R., Ma Y., Lyakhov A.O., Valle M., Gatti C. (2010)', ...
        'Evolutionary crystal structure prediction as a method', ...
        'for the discovery of minerals and materials.', ...
        'Rev. Mineral. Geochem. 71, 271-298' ...
        };
    formatted_rows = createHeader_wrap(text, 'left');
    for i=1:size(formatted_rows,2)
        fprintf(fp,[formatted_rows{i} '\n']);
    end
end

fprintf(fp,'            Job Starts at       %30s\n', datestr(now));
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for system description') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp,'                        Dimensionality  :  %2d\n', ORG_STRUC.dimension);
fprintf(fp,'                        Molecular       :  %2d (1:Yes, 0,No)\n', ORG_STRUC.molecule);
fprintf(fp,'                   Variable Composition :  %2d (1:Yes, 0,No)\n', ORG_STRUC.varcomp);

fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for atomic description') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);

fprintf(fp,'    There are %1d types of atoms in the system: ', length(ORG_STRUC.atomType));
for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'%5s', megaDoof(ceil(ORG_STRUC.atomType(i))));
end

fprintf(fp,'\n');

for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'    Minimum distances:                 %5s: ', megaDoof(ceil(ORG_STRUC.atomType(i))));
    for j=1:length(ORG_STRUC.atomType)
        fprintf(fp,'%4.2f  ', ORG_STRUC.minDistMatrice(i,j));
    end
    fprintf(fp,'\n');
end
fprintf(fp,'\n');
for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'           Good Bonds:                 %5s: ', megaDoof(ceil(ORG_STRUC.atomType(i))));
    for j=1:length(ORG_STRUC.atomType)
        fprintf(fp,'%4.2f  ', ORG_STRUC.goodBonds(i,j));
    end
    fprintf(fp,'\n');
end
fprintf(fp,'\n');
fprintf(fp,'            valences                        : ');
for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'%4.2f  ', ORG_STRUC.valences(i));
end
fprintf(fp,'\n');


if ORG_STRUC.molecule==1
    fprintf(fp,'    There are%2d types of molecules in the system: \n', size(ORG_STRUC.numMols,2));
    for i=1:size(ORG_STRUC.numIons,2)
        STDMOL = ORG_STRUC.STDMOL(i);
        format = STDMOL.format;
        Zmatrix = STDMOL.ZMATRIX;
        fprintf(fp,'%6s \n', STDMOL.name);
        fprintf(fp,'Please see the MOL_* files for the details.\n');
        fprintf(fp,'The calculated Zmatrix is:\n');
        fprintf(fp,'Atom Bond-length Bond-angle Torsion-angle  i  j  k\n');
        fprintf(fp,'      (Angstrom)  (Degree)    (Degree)   \n');
        for j = 1:length(STDMOL.types)
            type = megaDoof(ceil(ORG_STRUC.atomType(STDMOL.types(j))));
            Zmatrix(j,2) = Zmatrix(j,2)*180/pi;
            Zmatrix(j,3) = Zmatrix(j,3)*180/pi;
            fprintf(fp,'%2s   %10.4f %10.4f %10.4f   %4d %2d %2d\n', type, Zmatrix(j,:), format(j,:));
        end
    end
    fprintf(fp,'\n');
end

fprintf(fp,'    The investigated system is: ');
if ORG_STRUC.varcomp ==1
    for i=1:size(ORG_STRUC.numIons,1)
        fprintf(fp,'[');
        for j=1:size(ORG_STRUC.numIons,2)
            if ORG_STRUC.numIons(i,j)>0
                if ORG_STRUC.molecule==1
                    fprintf(fp,'%6s_%1d', ORG_STRUC.STDMOL(j).name, ORG_STRUC.numIons(i,j));
                else
                    fprintf(fp,'%2s_%1d', megaDoof(ceil(ORG_STRUC.atomType(j))), ORG_STRUC.numIons(i,j));
                end
            end
        end
        fprintf(fp,']');
        if i<size(ORG_STRUC.numIons,1)
            fprintf(fp,' --- ');
        end
    end
elseif ORG_STRUC.molecule==1
    for i=1:size(ORG_STRUC.numMols,2)
        fprintf(fp,'%6s_%1d', ORG_STRUC.STDMOL(i).name, ORG_STRUC.numMols(i));
    end
else
    for i=1:length(ORG_STRUC.atomType)
        fprintf(fp,'%2s_%2d  ', megaDoof(ceil(ORG_STRUC.atomType(i))), ORG_STRUC.numIons(i));
    end
end
fprintf(fp,'\n');

if ORG_STRUC.dimension==2
    fprintf(fp,'The reconstruction cell can be varied up to %2d multiplications: \n', ORG_STRUC.reconstruct);
    fprintf(fp,'Please see the POSCAR_SUBSTRATE for the details of substrate environment.\n');
    fprintf(fp,'\n');
end

fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for evolutionary algorithm') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp,'                 Number of Generations  :    %4d\n', ORG_STRUC.numGenerations);
fprintf(fp,'               Initial Population Size  :    %4d\n', ORG_STRUC.initialPopSize);
fprintf(fp,'               General Population Size  :    %4d\n', ORG_STRUC.populationSize);


fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('AB INITIO CALCULATIONS') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
if ~isempty(ORG_STRUC.ExternalPressure)
    fprintf(fp,'*  External Pressure is : %6.4f  GPa*\n', ORG_STRUC.ExternalPressure);
end
if ~ORG_STRUC.constLattice
    if ORG_STRUC.varcomp == 1   %for varcomp, we need to output volume for each block
        if ORG_STRUC.molecule == 1 %To be reconsidered later: for the case of varcomp of CaCl2 + H2O (2 blocks, 3 MOLs)
            %uff = ORG_STRUC.STDMOL;
            %fprintf(fp, '\n* Initial Volume for MOls : \n' );
            %for i = 1:length(uff)
            %	fprintf(fp, '     Volume of MOl_%d : %6.3f  A^3\n', i, ORG_STRUC.latVolume(i));
            %end
        else
            fprintf(fp, '\n*  Initial Volume for atom/block : \n' );
            for i = 1:size(ORG_STRUC.numIons,1)
                fprintf(fp, '    Volume of atom/block  ');
                for j = 1:length(ORG_STRUC.numIons(i,:) )
                    if ORG_STRUC.numIons(i,j)~=0
                        fprintf(fp, ' %s_%d', megaDoof(ORG_STRUC.atomType(j)), ORG_STRUC.numIons(i,j));
                    end
                end
                fprintf(fp, '  : %6.3f   A^3\n', ORG_STRUC.latVolume(i) );
            end
        end
    else %for fixcomp, we need to output the single volume value
        fprintf(fp, '     Estimated Volume : %6.3f  A^3\n', ORG_STRUC.latVolume);
    end
elseif ORG_STRUC.dimension == 3  % bulk + fixed_lattice
    fprintf(fp, '\n* This is a fixed lattice calculation \n' );
    for i = 1:3
        fprintf(fp, '%6.3f  %6.3f  %6.3f \n', ORG_STRUC.lattice(i,:) );
    end
end
fprintf(fp, '\n');

total_step = length(ORG_STRUC.abinitioCode);
fprintf(fp,'*  There are %2d local relaxation steps for each individual structure  *\n', total_step);
if ORG_STRUC.platform==0
    fprintf(fp,'Step  Abinitio Code                  Execute Command            K-resolution \n');
else
    fprintf(fp,'Step  Abinitio Code    K-resolution \n');
end
for step=1:total_step
    i = ORG_STRUC.abinitioCode(step);
    if i==0
        code=' no relaxation';
    elseif i==1
        code='     VASP     ';
    elseif i==2
        code='    SIESTA    ';
    elseif i==3
        code='     GULP     ';
    elseif i==4
        code='    LAMMPS    ';
    elseif i==5
        code=' Neu Network  ';
    elseif i==6
        code='   DMACRYS    ';
    elseif i==7
        code='    CP2K      ';
    elseif i==8
        code='    PWSCF     ';
    elseif i==9
        code='    FHI_aims  ';
    elseif i==10
        code='     ATK      ';
    elseif i==11
        code='   CASTEP     ';
    elseif i==12
        code='   Tinker     ';
    elseif i==13
        code='    MOPAC     ';
    elseif i==14
        code='  BoltzTraP   ';
    else
        fprintf(fp,'********************************************************************\n');
        fprintf(fp,'**      ERROR:    The code you selected is not vaild!             **\n');
        fprintf(fp,'**                 Program STOPS!!!!!!!!!!!                       **\n');
        fprintf(fp,'********************************************************************\n');
    end
    if sum(i==[3 4 5 6])==1 %for GULP, LAMMPS, NN, and DMACRYS
        ORG_STRUC.Kresol(step)=0;
    elseif isempty(ORG_STRUC.Kresol(step))
        fprintf(fp,'********************************************************************\n');
        fprintf(fp,'**      ERROR:   You did not specify Kresolution for each Step    **\n');
        fprintf(fp,'**                 Program STOPS!!!!!!!!!!!                       **\n');
        fprintf(fp,'********************************************************************\n');
    end
    if ORG_STRUC.platform==0
        fprintf(fp,' %2d %12s  %36s  %12.3f\n', step, code, ORG_STRUC.commandExecutable{step}, ORG_STRUC.Kresol(step));
    else
        fprintf(fp,' %2d %12s  %12.3f\n', step, code, ORG_STRUC.Kresol(step));
    end
end
fprintf(fp,'\n');

if ORG_STRUC.platform==0
    fprintf(fp,'The calculations are performed in nonParallel mode on the local machine\n');
elseif ORG_STRUC.platform==1
    fprintf(fp,'The script for job submission is prepared seperately in Submission/*_local.m\n');
elseif ORG_STRUC.platform==2
    fprintf(fp,'The script for job submission is prepared seperately in Submission/*_remote.m\n');
elseif ORG_STRUC.platform==3
    fprintf(fp,'The calculations are performed in Parallel mode on CFN supercomputer\n');
elseif ORG_STRUC.platform==4 | ORG_STRUC.platform==5
    fprintf(fp,'The calculations are performed in Parallel mode on QSH supercomputer\n');
elseif ORG_STRUC.platform==6
    fprintf(fp,'The calculations are performed in Parallel mode on xservDE supercomputer\n');
elseif ORG_STRUC.platform==7
    fprintf(fp,'The calculations are performed in Parallel mode on MIPT supercomputer\n');
end
fprintf(fp,'%4d parallel calculations are performed simutaneously\n', ORG_STRUC.numParallelCalcs);
fprintf(fp,'\n');

fclose(fp);
