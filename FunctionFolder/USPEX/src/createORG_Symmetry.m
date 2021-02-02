function createORG_Symmetry(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%[nothing, doSpaceGroup] = unix (['./getStuff ' inputFile ' doSpaceGroup 1']);
doSpaceGroup = python_uspex(getPy, ['-f ' inputFile ' -b doSpaceGroup -c 1']);
if isempty(doSpaceGroup)
    if (ORG_STRUC.dimension==0) | (ORG_STRUC.dimension==2) | (ORG_STRUC.dimension==-3)
        doSpaceGroup = '0'; % default
    else
        doSpaceGroup = '1'; % default
    end
end
ORG_STRUC.doSpaceGroup = str2num(doSpaceGroup);


% do we want to symmetrize the structure using Stokes SG determination code (using symmetrized.cif)
% done at the last optimization step (when not shaking - symmetrize)
%[nothing, symmetrize] = unix(['./getStuff ' inputFile ' symmetrize 1']);
symmetrize = python_uspex(getPy, ['-f ' inputFile ' -b symmetrize -c 1']);
if ~isempty(symmetrize)
    ORG_STRUC.symmetrize = str2num(symmetrize);
end
% Space group determination tolerance
%[noathing, SGtolerance] = unix(['./getStuff ' inputFile ' SymTolerance 1']);
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
    SGtolerance = deblank(SGtolerance);
    if strcmp(lower(SGtolerance), 'high')
        ORG_STRUC.SGtolerance = 0.04;
    elseif strcmp(lower(SGtolerance), 'medium')
        ORG_STRUC.SGtolerance = 0.08;
    elseif strcmp(lower(SGtolerance), 'low')
        ORG_STRUC.SGtolerance = 0.15;
    else % number!
        ORG_STRUC.SGtolerance = str2num(SGtolerance);
    end
end
if isempty(ORG_STRUC.SGtolerance) %In case some illegal strings
    ORG_STRUC.SGtolerance = 0.08;
end
% what point symmetries have to be satisfied by clusters created randomly
% for crystals - what symmetry groups the randomly created crystals should belong to, format : 1-3 5 6-9 15  etc
%[nothing, nsym] = unix(['./getFromTo symmetries endSymmetries ' inputFile]);
nsym = python_uspex(getPy, ['-f ' inputFile ' -b symmetries -e EndSymmetries']);
if isempty(nsym)
    if ORG_STRUC.dimension==0
        nsym = 'E C2 D2 C4 C3 C6 T S2 Ch1 Cv2 S4 S6 Ch3 Th Ch2 Dh2 Ch4 D3 Ch6 O D4 Cv3 D6 Td Cv4 Dd3 Cv6 Oh Dd2 Dh3 Dh4 Dh6 Oh C5 S5 S10 Cv5 Ch5 D5 Dd5 Dh5 I Ih '; % default
    elseif ORG_STRUC.dimension==3
        nsym = '2-230'; % default
    elseif ORG_STRUC.dimension==-2 | ORG_STRUC.dimension==1 | ORG_STRUC.dimension==-3
        nsym = '2-17';
    end
end
if ORG_STRUC.dimension==0
    nsym(end) = [];
    ORG_STRUC.nsym = nsym;
    c1 = findstr(nsym, ' ');
    c = sort(str2num(['0 ' num2str(c1)]));
    c(end+1) = length(nsym) + 1;
    ind = zeros(1,2);
    indN = 0;
    for i = 2 : length(c)
        if c(i-1)+1 > c(i)-1
            continue
        end
        indN = indN + 1;
        ind(1,1) = c(i-1)+1;
        ind(1,2) = c(i)-1;
        if indN == 1
            ORG_STRUC.nsymN = ind;
        else
            ORG_STRUC.nsymN = cat(1,ORG_STRUC.nsymN,ind);
        end
    end
    if isempty(ORG_STRUC.nsym)
        ORG_STRUC.nsym = 'E';
        ORG_STRUC.nsymN = [1 1];
    end
else
    ORG_STRUC.nsym = zeros(1,230);
    c1 = findstr(nsym, ' ');
    c2 = findstr(nsym, '-');
    c = sort(str2num(['0 ' num2str(c1) ' ' num2str(c2)]));
    c(end+1) = length(nsym) + 1;
    ind1 = 1;
    for i = 2 : length(c)
        if c(i-1)+1 > c(i)-1
            continue
        end
        ind2 = str2num(nsym(c(i-1)+1 : c(i)-1));
        if ind2 == 0
            ind1 = 1;
            continue
        end
        if ~isempty(find(c2 == c(i-1)))
            for j = ind1 : ind2
                ORG_STRUC.nsym(j) = 1;
            end
        else
            ORG_STRUC.nsym(ind2) = 1;
        end
        ind1 = ind2;
    end
    if sum(ORG_STRUC.nsym) == 0
        ORG_STRUC.nsym(1) = 1;
    end
    ORG_STRUC.nsymN = [0 0];
end
% coefficient between mindist and symmetrization distance (1 by default, sometimes > 1 needed)
%[nothing, sym_coef] = unix(['./getStuff ' inputFile ' constraint_enhancement 1']);
sym_coef = python_uspex(getPy, ['-f ' inputFile ' -b constraint_enhancement -c 1']);
if ~isempty(sym_coef)
    ORG_STRUC.sym_coef = str2num(sym_coef);
end

