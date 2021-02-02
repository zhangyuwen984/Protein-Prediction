%This function write files that are used in the property calculation : 
% BoltzTraP calculation requires:
% CASE.struct   : structure info
% CASE.energy   : energy bands
% CASE.intrans  : input parameters
% BoltzTraP.def : format file
% Dong Dong 

% function writeProperty()
function Write_BoltzTraP(Ind_No)

    global ORG_STRUC;
%%%% STEP 1 preparation generic input file for boltztrap

% if additional programs in Property folder needed for property calculation 
% then copy them to the calculation folder
% [nothing, nothing] = unix('cp ../Property/w_struct .');
% [nothing, nothing] = unix('cp ../Property/BoltzTraP .');

% detemine the current folder name
  currentdir = pwd;
  [upperpath,dir] = fileparts(currentdir);


% USPEX rename the VASP results, it is necessary to recover the results from the last step
[nothing, nothing] = unix('cp CONTCAR_old CONTCAR');
[nothing, nothing] = unix('cp OUTCAR_old OUTCAR');
[nothing, nothing] = unix('cp OSZICAR_old OSZICAR');
[nothing, nothing] = unix('cp POSCAR_old POSCAR');

% delete previous results
  delete([dir '-*']);  %% DANGEROUS
  
  % Commented out Dong Dong's implementation. Fei Qi, 2015/07/27
  % Using the vasp2boltz.py script
  [nothing, USPEXPATH]=unix('echo -e $USPEXPATH');
  if(length(USPEXPATH) <= 1)
      USPEXPATH = ORG_STRUC.homePath
  end
  python_uspex([USPEXPATH(1:end-1) '/FunctionFolder/Tool/vasp2boltz.py'], num2str(Ind_No));
  
% % write CalcFoldX.struct
% % fortran code 'w_struct' takes care all
% [nothing, nothing] = unix(['$USPEXPATH/FunctionFolder/Tool/w_struct > ' dir '.struct']);


% % write CalcFoldX.energy
%  fp = fopen('EIGENVAL');

% % read fermi energy in eV
%  [tmp,E_fermi] = unix('./getStuff OUTCAR E-fermi 4');
%  E_fermi = str2num(E_fermi);

% % get the 6th lines
%  for i = 1 : 6
%    tmp = fgetl(fp);
%  end
%  line6 = str2num(tmp);

%  numEln = line6(1); % number of electrons
%  numKpt = line6(2); % number of k-points
%  numBnd = line6(3); % number of bands

% % if EIGENVAL contains no data, do not create case.energy file
% if ~feof(fp)
    
%     fpo = fopen([dir '.energy'],'w');

%     %line 1: comment
%     fprintf(fpo,'1.2 BoltzTraP \n');

%     %line 2: number of K-points
%     fprintf(fpo,'%d \n',numKpt);

%     %lines from 3rd
%     for i = 1 : numKpt
%         Kpt = fscanf(fp,'%g%g%g%g',[4]);
%         % k vector
%         fprintf(fpo,'  %g   %g   %g   %d\n', Kpt(1:3),numBnd);
%         % energy levels at k-vector
%         for j = 1 : numBnd
%             En=fscanf(fp,'%g%g',[2]);
%             % boltztrap requires Ry=13.605698eV; move Efermi to 0
%             fprintf(fpo,'   %10.9g \n',(En(2) - E_fermi)/13.605698);
%         end
%     end

%     fclose(fpo);
% end

%  fclose(fp); 

% % step 1.3
% % write CalcFoldX.intrans

%  fpo = fopen([dir '.intrans'],'w');
%  fprintf(fpo,'GENE                          # use generic interface\n');
%  fprintf(fpo,'0      0      0   0.0         # iskip (not presently used) idebug setgap shiftgap \n');
%  fprintf(fpo,'0.0000 0.0005 0.4   %6.1f     # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n', numEln);
%  fprintf(fpo,'CALC                          # CALC (calculate expansion coeff), NOCALC read from file\n');
%  fprintf(fpo,'5                             # lpfac, number of latt-points per k-point\n');
%  fprintf(fpo,'BOLTZ                         # run mode\n');
%  fprintf(fpo,'.15                           # (efcut) energy range of chemical potential\n');
%  fprintf(fpo,'800.    50.                   # Tmax, temperature grid\n');
%  fprintf(fpo,'-1.                           # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)\n');
%  fprintf(fpo,'HISTO\n');
%  fprintf(fpo,'0 0 0 0 0\n');
%  fprintf(fpo,'2\n');
%  fprintf(fpo,'1E20 -1E20\n');
%  fclose(fpo);


% % write a BoltzTraP file dictionary
%  fpo=fopen('BoltzTraP.def', 'w');

%  fprintf(fpo,['5,  ' dir '.intrans,      ''old'',      ''formatted'' ,0 \n']);
%  fprintf(fpo,['6,  ' dir '.outputtrans,  ''unknown'',  ''formatted'' ,0 \n']);
%  fprintf(fpo,['20, ' dir '.struct,       ''old'',      ''formatted'' ,0 \n']);
%  fprintf(fpo,['10, ' dir '.energy,       ''old'',      ''formatted'' ,0 \n']);
%  fprintf(fpo,['48, ' dir '.engre,        ''unknown'',  ''unformatted'',0 \n']);
%  fprintf(fpo,['49, ' dir '.transdos,     ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['50, ' dir '.sigxx,        ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['51, ' dir '.sigxxx,       ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['21, ' dir '.trace,        ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['22, ' dir '.condtens,     ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['24, ' dir '.halltens,     ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['30, ' dir '_BZ.cube,      ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['35, ' dir '_band.dat,     ''unknown'',  ''formatted'',0 \n']);
%  fprintf(fpo,['36, ' dir '_band.gpl,     ''unknown'',  ''formatted'',0 \n']);
 
%  fclose(fpo);

