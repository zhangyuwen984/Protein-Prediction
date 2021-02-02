function jobNumber = submitJob_remote(USPEX, Index, Ind_No)

global POP_STRUC
global ORG_STRUC

%-------------------------------------------------------------
%This routine is to check if the submitted job is done or not
%2   : whichCluster (default 0, 1: local submission; 2: remote submission)
%C-20GPa : remoteFolder
%-------------------------------------------------------------

%-------------------------------------------------------------
%Step1: To prepare the job script, runvasp.sh
  fp = fopen('myrun', 'w');
  fprintf(fp, '#!/bin/bash\n');
  fprintf(fp, '#SBATCH -N 1\n');
  fprintf(fp, '#SBATCH -n 1\n');
  fprintf(fp, '#SBATCH -p RT\n');
  fprintf(fp, '#SBATCH --job-name=1shf\n');
  fprintf(fp, '#SBATCH -o out -e err\n');
  fprintf(fp, '#SBATCH --comment="calculate protein ishf"\n');
  fclose(fp);
%------------------------------------------------------------------------------------------------------------
%Step 2: Copy the files to the remote machine

%Step2-1: Specify the PATH to put your calculation folder
Home = ['/home/AI/shipilov.ab/statistic_collection/1cei']; %'pwd' of your home directory of your remote machine
Address = 'shipilov.ab@calc.cod.phystech.edu'; %your target server: username@address
Path = [Home '/' USPEX '/CalcFold' num2str(Index)];  %Just keep it

%------------------------------------------------------------------------------------------------------------
%Step 3: to submit the job and get JobID, i.e. the exact command to submit job.

%Здесь твой код, который ты хочешь запустить в каждом CalcFold в формате: 
[nothing, nothing] = unix(['cat /dev/null > myrun']);
[nothing, nothing] = unix(['echo "#!/bin/sh"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -o out"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -p RT"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -J U-' num2str(POP_STRUC.generation),'I',num2str(Ind_No),'S',num2str(POP_STRUC.POPULATION(Ind_No).Step),'"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -t 06:00:00"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -N 1"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH -n 1"  >> myrun']);
[nothing, nothing] = unix(['echo "#SBATCH --comment=\"1cei protein calculation\""  >> myrun']);
[nothing, nothing] = unix(['echo cd ' Path ' >> myrun']);
%[nothing, nothing] = unix(['echo "mpirun vasp_std> log" >> myrun']);
[a,b]=unix(['echo "~/tools/miniconda3/envs/env/bin/python ' Path '/random_protein.py input 1 pseudo memory > output"  >> myrun']);

%Создаем папку USPEX в директории Home, а в ней папки CalcFold
try
[a,b]=unix(['ssh ' Address ' "cd ' Home '; mkdir ' USPEX '"' ]);  
catch
end

try
[a,b]=unix(['ssh ' Address ' "cd ' Home '; cd ' USPEX '; mkdir ' Path '"' ]);
catch
end

%Копируем все необходимые файлы из папки на рюрике в папку на миптовском кластере, в которой будет исполняться код
[nothing, nothing] = unix(['scp -r * ' Address ':' Path]);

[a,v]=unix(['ssh ' Address ' "cd ' Path '; sbatch myrun"']);
start_marker=findstr(v,'job ');
jobNumber = v(start_marker(1)+4:end-1);
disp([ 'Individual : ' num2str(Ind_No) ' -- JobID :', num2str(jobNumber) ]);
