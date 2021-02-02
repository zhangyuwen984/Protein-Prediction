function createORG_Surface(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%for surfaces

%[nothing, reconstruct] = unix (['./getStuff ' inputFile ' reconstruct 1']);
reconstruct = python_uspex(getPy, ['-f ' inputFile ' -b reconstruct -c 1']);
if ~isempty(reconstruct)
    ORG_STRUC.reconstruct = str2num(reconstruct);
end
%[nothing, thicknessS] = unix (['./getStuff ' inputFile ' thicknessS 1']);
%thicknessS = python_uspex(getPy, ['-f ' inputFile ' -b thicknessS -c 1']);
%if ~isempty(thicknessS)
%    ORG_STRUC.thicknessS = str2num(thicknessS);
%end

%[nothing, thicknessB] = unix (['./getStuff ' inputFile ' thicknessB 1']);
%thicknessB = python_uspex(getPy, ['-f ' inputFile ' -b thicknessB -c 1']);
%if ~isempty(thicknessB)
%    ORG_STRUC.thicknessB = str2num(thicknessB);
%end

%[nothing, bulk_stoi] = unix(['./getFromTo StoichiometryStart StoichiometryEnd ' inputFile]);
bulk_stoi = python_uspex(getPy, ['-f ' inputFile ' -b StoichiometryStart -e StoichiometryEnd']);
if isempty(bulk_stoi)
    ORG_STRUC.bulk_stoi = [1 1];
else
    ORG_STRUC.bulk_stoi = str2num(bulk_stoi);
end

%[nothing, E_AB] = unix (['./getStuff ' inputFile ' E_AB 1']);
E_AB = python_uspex(getPy, ['-f ' inputFile ' -b E_AB -c 1']);
if ~isempty(E_AB)
    ORG_STRUC.E_AB = str2num(E_AB);
end

%[nothing, E_A] = unix (['./getStuff ' inputFile ' Mu_A 1']);
E_A = python_uspex(getPy, ['-f ' inputFile ' -b Mu_A -c 1']);
if ~isempty(E_A)
    ORG_STRUC.E_A = str2num(E_A);
end

%[nothing, E_B] = unix (['./getStuff ' inputFile ' Mu_B 1']);
E_B = python_uspex(getPy, ['-f ' inputFile ' -b Mu_B -c 1']);
if ~isempty(E_B)
    ORG_STRUC.E_B = str2num(E_B);
end


