function PickUp()

global ORG_STRUC
global POP_STRUC
global POOL_STRUC
global USPEX_STRUC

if ORG_STRUC.pickUpFolder == 0
    ORG_STRUC.pickUpFolder = str2num(ORG_STRUC.resFolder(8:end))-1;
end

cd (['results' num2str(ORG_STRUC.pickUpFolder)])
try
    load('ANTISEEDS.mat')
catch
end

if ORG_STRUC.pickUpGen == 0
    ORG_STRUC.pickUpGen = 1;
    while exist (['generation' num2str( ORG_STRUC.pickUpGen)]) == 7
        ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen + 1;
    end
    ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen - 1;
end

cd (['generation' num2str(ORG_STRUC.pickUpGen)])
[nothing, nothing] = unix('pwd');
load('POP_STRUC.mat')
load('USPEX.mat')
if exist('POOL.mat')
    load('POOL.mat')
end

cd ../

%Поиск структуры с самым большим номером среди того поколения, с которого стартуем
allStructs=[POP_STRUC.POPULATION().Number];
maxValue = max(allStructs);

%Удаление лишних структур из gatheredPDB, gatheredMAKE (тех, что относятся к поколениям после pickUpGen)
% и копирование его в новую папку
[nothing, nothing] = unix(['cp gatheredPDB >> gatheredPDB_copy']);
%[nothing, nothing] = unix(['cat all_gatheredPDB >> tmp_gatheredPDB']);
[nothing, nothing] = unix(['cat old_gatheredPDB >> tmp_gatheredPDB']);
[nothing, nothing] = unix(['mv tmp_gatheredPDB  gatheredPDB_old']);
[nothing, nothing] = unix(['remove_excess.py gatheredPDB EA' num2str(maxValue + 1)]);
[nothing, nothing] = unix(['mv gatheredPDB  ../' num2str(ORG_STRUC.resFolder)]);

[nothing, nothing] = unix(['mv gatheredMAKE  gatheredMAKE_old']);
[nothing, nothing] = unix(['remove_excess.py gatheredMAKE EA' num2str(maxValue + 1)]);
[nothing, nothing] = unix(['mv gatheredMAKE  ../' num2str(ORG_STRUC.resFolder)]);

cd (['generation' num2str(ORG_STRUC.pickUpGen)])

if ORG_STRUC.fixRndSeed > 0
    rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
end

cd ../..


ORG_STRUC.initialPopSize = length(POP_STRUC.POPULATION);
ORG_STRUC.pickUpNCount = POP_STRUC.bodyCount;

POP_STRUC.resFolder = ORG_STRUC.resFolder;
ORG_STRUC.pickedUP = 1;
disp('This Calculations has been picked up from an older calculation');
safesave ([ORG_STRUC.resFolder '/PickedUP_POP_STRUC.mat'],POP_STRUC)

ORG_STRUC.numGenerations = ORG_STRUC.numGenerations + POP_STRUC.generation;

safesave([ORG_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
try
    safesave([ORG_STRUC.resFolder '/POOL.mat'], POOL_STRUC);
catch
end

