function createORGStruc(inputFile)

% USPEX Version 9.4.1
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
warning off
global ORG_STRUC

ORG_STRUC = struct('tournament',{},'dimension',{},'varcomp',{},'molecule',{},...%EA
    'PhaseDiagram',{}, 'maxErrors',{},'numProcessors', {}, 'numParallelCalcs', {},... %others
    'pluginType',{}, 'currentGenDone',{}, 'startNextGen',{},...                             %Pop and Generation
    'stopCrit',{}, 'reoptOld', {}, 'bestFrac',{},'repeatForStatistics',{},...
    'fracGene',{},'fracRand',{},'fracPerm',{}, 'fracLatMut', {},...        %fraction
    'fracRotMut',{},'fracAtomsMut',{}, 'AutoFrac',{}, 'stopFitness', {},...
    'howManyAtomMutations',{}, 'howManyRand',{},...                         %how many
    'topologyRandom',{},'coordinatingNumbers',{},...                                %topology
    'howManyOffsprings',{},'maxDistHeredity',{},...             %heredity
    'minSlice', {}, 'maxSlice', {}, 'percSliceShift', {},...
    'softMutOnly', {}, 'softMutTill', {},...                             %softmutation
    'chargeNeutrality',{}, 'checkConnectivity',{}, 'howManyMut',{},...
    'valences',{}, 'NvalElectrons', {}, 'goodBonds',{}, 'mutationRate',{},...
    'howManyPermutations', {},'howManySwaps', {}, 'specificSwaps', {},...  %permutation
    'fracSecSwitch',  {}, 'howManySecSwitch',   {}, ... % Secondary Switch operator for proteins
    'fracShiftBorder',{}, 'howManyShiftBorder', {}, ... % Shift Border operator for proteins
    'abinitioCode',{}, 'Kresol',{}, 'commandExecutable',{},...   %ab-initio calculation
    'remoteFolder', {}, 'platform',{}, 'ExternalPressure', {},...
    'latVolume',{},'lattice',{},'constLattice',{},'splitInto',{},...  %lat and vol
    'homePath', {}, 'USPEXPath',{}, 'specificFolder', {}, 'log_file', {},     ... %path and output
    'atomType',{}, 'numIons', {},                        ...                %atom
    'averageFitness',{},'correlation_coefficient',{},... %aging
    'antiSeedsMax',{}, 'antiSeedsSigma',{}, 'ordering',{}, 'antiSeedsActivation',{},...
    'optType',{}, 'opt_sign',{},'cor_dir', {}, ...                        %optimization
    'RmaxFing',{}, 'deltaFing', {}, 'sigmaFing', {}, 'doFing', {},...     %Fingerprint
    'toleranceFing',{}, 'toleranceBestHM',{}, 'dynamicalBestHM',{}, ...
    'keepBestHM',{}, 'alpha',{}, 'erf_table',{}, 'weight',{},...
    'numMols',{},'STDMOL',{}, 'checkMolecules', {},...                %Molecule
    'thicknessS',{}, 'thicknessB',{},'bulk_lat',{},'bulk_pos',{},...      %Surface
    'bulk_ntyp',{},'bulk_atyp',{},'reconstruct',{},...
    'bulk_stoi',{}, 'E_AB',{},'E_A',{},'E_B',{},...
    'vacuumSize',{}, 'numParents',{}, 'manyParents',{}, ...              %Nanoparticles
    'firstGeneSplit',{},'firstGeneSplitAll',{},'splitN',{},'minAt',{},'maxAt',{},'firstGeneMax',{},... %varcomp;
    'howManyTrans',{}, 'fracTrans',{}, 'specificTrans',{}, ...           %transmuation
    'doSpaceGroup',{}, 'nsym',{},'nsymN',{},...                              %Symmetry
    'sym_coef',{}, 'symmetrize',{}, 'SGtolerance',{},...                     %Symmetry
    'pickUpYN',{},'pickUpGen',{}, 'pickUpFolder',{},'pickUpNCount',{},...     %Pick up
    'minAngle',{}, 'minDiagAngle',{}, 'minDistMatrice',{},...              %Constraints
    'CenterminDistMatrice',{}, 'minVectorLength',{},...                    %Constraints
    'numGenerations',{}, 'populationSize',{}, 'initialPopSize',{},...      %Pop and Generation
    'spin',{}, 'ldaU',{}, 'magRatio',{}, 'fracSpin',{}, 'howManySpinmutations',{}, ...      %for-SPIN-function
    'fixRndSeed',{}, 'collectForces',{});

% defaults
createORGDefault();
% %%%%%%%% RESULTS
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
    folderNum = folderNum + 1;
    ORG_STRUC.resFolder = ['results' num2str(folderNum)];
    [nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'w');


% Create a header automatically:
description    = 'Evolutionary Algorithm Code for Structure Prediction';
funcFold       = [ORG_STRUC.USPEXPath '/FunctionFolder'];
formatted_rows = createHeader('USPEX_pic.txt', description, funcFold);

for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end


% Cite:
text = {'Please cite the following suggested papers',       ...
    'when you publish the results obtained from USPEX:' ...
    };
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {'Oganov A.R., Glass C.W. (2006). Crystal structure prediction', ...
    'using evolutionary algorithms: principles and applications.',  ...
    'J. Chem. Phys. 124, 244704', ...
    '', ...
    'Oganov A.R., Stokes H., Valle M. (2011)', ...
    'How evolutionary crystal structure prediction works - and why.', ...
    'Acc. Chem. Res. 44, 227-237', ...
    '', ...
    'Lyakhov A.O., Oganov A.R., Stokes H., Zhu Q. (2013)', ...
    'New developments in evolutionary structure prediction algorithm USPEX.', ...
    'Comp. Phys. Comm., 184, 1172-1182' ...
    };
formatted_rows = createHeader_wrap(text, 'left');
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {datestr(now)};
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end
fprintf(fp,'\n');


%%%%%%% GET the code and cluster.
createORG_System(inputFile);
createORG_AbinitCode(inputFile);
createORG_Symmetry(inputFile);
createORG_EA(inputFile);
if ORG_STRUC.dimension == 2 | ORG_STRUC.dimension == -3
    if exist('POSCAR_SUBSTRATE')
        read_SUBSTRATE();
        createORG_Surface(inputFile);
    else
        disp('Please provide the POSCAR_SUBSTATE file');
        disp('Otherwise, the calculation can not start.');
        disp('The programs STOPS here!');
        fprintf(fp, 'Please provide the POSCAR_SUBSTATE file.\n');
        fprintf(fp, 'Otherwise, the calculation can not start.\n');
        fprintf(fp, 'The programs STOPS here!\n');
        quit
    end
    
end

if exist('MOL_1')
    createORG_Molecule(inputFile);
end

if (ORG_STRUC.molecule == 1) & ~exist('MOL_1')
    disp(['For Molecules, you did not give MOL_x files']);
    disp(['The program is STOPS here........']);
    fprintf(fp, 'For Molecules, you did not give MOL_x files.\n');
    fprintf(fp, 'The program STOPS here......... \n');
    quit
end
fclose(fp);
createORG_Fingerprint(inputFile);

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LATTICE manipulations %%%%%%%%%%%% MOVE OUT FROM SYSTEM
% varcomp does not support fixed lattice
% information on the lattice
%[nothing, lattice] = unix(['./getFromTo Latticevalues Endvalues ' inputFile]);
%ORG_STRUC.lattice = str2num(lattice);
lattice_INPUT = python_uspex(getPy, ['-f ' inputFile ' -b Latticevalues -e Endvalues'], 1);

if size(lattice_INPUT, 1) == 1
    if (size(lattice_INPUT, 2) == 6)     % fixed lattice in a form {a b c alpha beta gamma}?
        lattice_INPUT(4:6) = lattice_INPUT(4:6)*pi/180;         % Matlab works with radians, input - in degrees
        ORG_STRUC.constLattice = 1;
        ORG_STRUC.lattice = latConverter(lattice_INPUT);   % 3x3 form
        ORG_STRUC.latVolume = det(ORG_STRUC.lattice);
    else  % varcomp
        ORG_STRUC.latVolume = lattice_INPUT;
        ORG_STRUC.constLattice = 0;
    end
elseif (sum(size(lattice_INPUT) == [3 3]) == 2)  %3*3 matrix
    ORG_STRUC.latVolume = det(lattice_INPUT);
    ORG_STRUC.constLattice = 1;
    ORG_STRUC.lattice = lattice_INPUT;
else %empty
    if ORG_STRUC.dimension~=2 & ORG_STRUC.dimension~=-3 & ORG_STRUC.dimension~=-4
        estimateInputVolume();
        ORG_STRUC.constLattice = 0;
    end
end

if ORG_STRUC.dimension==0 | ORG_STRUC.dimension==2 | ORG_STRUC.dimension==-3
    ORG_STRUC.constLattice = 1;
end
%%%%%%% END LATTICE manipulations %%%%%%%%%%%%
% how many individuals in the first generation
%[nothing, initialPopSize] =unix (['./getStuff ' inputFile ' initialPopSize 1']);
initialPopSize = python_uspex(getPy, ['-f ' inputFile ' -b initialPopSize -c 1']);
if isempty(initialPopSize)
    initialPopSize = num2str(ORG_STRUC.populationSize); % default
end
ORG_STRUC.initialPopSize = str2num(initialPopSize);

if ORG_STRUC.varcomp==1 & ORG_STRUC.dimension~=2
    % how many different compositions we calculate for first generation
    %[nothing, ans] = unix (['./getStuff ' inputFile ' firstGeneMax 1']);
    ans = python_uspex(getPy, ['-f ' inputFile ' -b firstGeneMax -c 1']);
    if isempty(ans)
        ORG_STRUC.firstGeneMax = round(ORG_STRUC.initialPopSize/4); % default
    else
        ORG_STRUC.firstGeneMax = str2num(ans);
    end
end

% number of atoms to split cell to
%[nothing, splitInto] = unix(['./getFromTo splitInto EndSplitInto ' inputFile]);
splitInto =  python_uspex(getPy, ['-f ' inputFile ' -b splitInto -e EndSplitInto']);
if ~isempty(splitInto)
    ORG_STRUC.splitInto = str2num(splitInto);
end

% how many times to repeat the calculations with the same initial conditions
% this option allows to make some statistics rerunning USPEX a few times after calculations are done
%[nothing, repeatForStatistics] = unix(['./getStuff ' inputFile ' repeatForStatistics 1']);
repeatForStatistics = python_uspex(getPy, ['-f ' inputFile ' -b repeatForStatistics -c 1']);
if ~isempty(repeatForStatistics)
    ORG_STRUC.repeatForStatistics = str2num(repeatForStatistics);
end
repeatForStatistics = num2str(ORG_STRUC.repeatForStatistics);
if ORG_STRUC.repeatForStatistics > 1
    if ~exist('multiple_runs')
        [nothing, nothing] = unix(['echo "' repeatForStatistics ' : runs_left" > multiple_runs']);
    end
end

% connectivity check in hardness, distanceCheck etc
%[nothing, checkConnectivity] = unix(['./getStuff ' inputFile ' checkConnectivity 1']);
checkConnectivity = python_uspex(getPy, ['-f ' inputFile ' -b checkConnectivity -c 1']);
if isempty(checkConnectivity)
    if ORG_STRUC.optType == 3
        ORG_STRUC.checkConnectivity=1;  %default
    else
        ORG_STRUC.checkConnectivity=0;  %default
    end
else
    ORG_STRUC.checkConnectivity = str2num(checkConnectivity);
end

%[nothing, PhaseDiagram] = unix(['./getStuff ' inputFile ' PhaseDiagram 1']);
%ORG_STRUC.PhaseDiagram = str2num(PhaseDiagram);
PhaseDiagram = python_uspex(getPy, ['-f ' inputFile ' -b PhaseDiagram -c 1']);
if ~isempty(PhaseDiagram)
    ORG_STRUC.PhaseDiagram = str2num(PhaseDiagram);
else
    if ORG_STRUC.dimension ~= 3
        ORG_STRUC.PhaseDiagram = 0 ;
        %	else
        %		ORG_STRUC.PhaseDiagram = 0 ;
    end
end

stopFitness = python_uspex(getPy, ['-f ' inputFile ' -b stopFitness -c 1']);
if ~isempty(stopFitness)
    ORG_STRUC.stopFitness = str2num(stopFitness);
end

%[nothing, thicknessS] = unix (['./getStuff ' inputFile ' thicknessS 1']);
thicknessS = python_uspex(getPy, ['-f ' inputFile ' -b thicknessS -c 1']);
if ~isempty(thicknessS)
    ORG_STRUC.thicknessS = str2num(thicknessS);
end
%[nothing, thicknessB] = unix (['./getStuff ' inputFile ' thicknessB 1']);
thicknessB = python_uspex(getPy, ['-f ' inputFile ' -b thicknessB -c 1']);
if ~isempty(thicknessB)
    ORG_STRUC.thicknessB = str2num(thicknessB);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeORG()
%checkORG()
%[nothing,homePath] = unix('pwd');
%homePath(end) = [];
%checkPOTCARS(ORG_STRUC.atomType, ORG_STRUC.abinitioCode, ORG_STRUC.var_comp, homePath)
