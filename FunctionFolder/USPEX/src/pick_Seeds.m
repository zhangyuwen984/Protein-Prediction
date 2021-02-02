function loops_seed = pick_Seeds()

% USPEX Version 9.3
% Change: block varcomp support added
% Change: the seeds structure would be added into the end of the POPULATION
% Change: VASP-5.2 format
global POP_STRUC
global ORG_STRUC
warning off

%atomType_Seeds = zeros(1, length(ORG_STRUC.atomType));


%-------------------------------------
cd Seeds

disp(' ');
disp('Read Seeds ... ')
disp(' ');


% reading POSCARS
loops_seed = 0;
if exist([ 'POSCARS_', num2str(POP_STRUC.generation) ])
    if exist('POSCARS')
       [a,b]=unix([ 'cat POSCARS >> POSCARS_' num2str(POP_STRUC.generation) ]);
    end
    [a,b]=unix([ 'cp POSCARS_' num2str(POP_STRUC.generation) ' POSCARS' ]);
end

if exist('POSCARS')
    [fid,message] = fopen('POSCARS');
    while 1
        try
            
            tmp = fgetl(fid);      % system description
            if tmp==-1; break; end  % checking the end of the Seeds file.
            
            numMolsRead  = numMolReadSeeds( tmp );   % Read the MOL user setting.
            
            scale_factor = fgetl(fid); % 1 = numbers in angstrems
            optlattice   = fscanf(fid,'%g',[3,3]); % lattice vectors
            candidate_lat= optlattice'*str2num(scale_factor);
            
            tmp = fgetl(fid);  % we don't need this line
            
            %-------------------------%
            %--  To get atom Type   --%
            %-------------------------%
            atomType = fgetl(fid); % VASP5, elemental type
            
            if str2num(atomType)
                USPEXmessage(551, '', 0);
                break;
            end
            
            ntyp = fgetl(fid); % types of atoms aka atomic numbers
            ntyp = str2num(ntyp);
            N_type = length(ntyp);
            atomType1 = GetElement(N_type, atomType);            

            %-- check the consistent of AtomType
            isAtomTypeOK=1;
            if N_type > length(ORG_STRUC.atomType) %more atomTypes
                USPEXmessage(552, '', 0);
                isAtomTypeOK = 0;
            elseif sum(ismember(atomType1, ORG_STRUC.atomType))<N_type 
                USPEXmessage(553, '', 0);
                isAtomTypeOK = 0;
            end
            
            if  isAtomTypeOK == 0
                warningStr = ['Seeds : Elemental types are inconsistent when reading Seeds-' num2str(loops_seed+1) ];
                USPEXmessage(0, warningStr, 0);
                break;
            end

            %if isAtomTypeOK==0; continue; end
            %-------------------------%
            %----    End AtomType   --%
            %-------------------------%
            
            numIons = zeros(1,length(ORG_STRUC.atomType));
            atomType_Seeds = atomType1;
            
            for i=1:length(ORG_STRUC.atomType)
                for j=1:length(ntyp)
                    if atomType1(j)==ORG_STRUC.atomType(i)
                        numIons(i) = ntyp(j);
                    end
                end
            end

            if ORG_STRUC.molecule == 1
                comp=zeros(size(ORG_STRUC.numIons,2), length(ORG_STRUC.atomType));
                for i=1:size(ORG_STRUC.numIons,2)
                    for j=1:length(ORG_STRUC.STDMOL(i).types)
                        for k=1:length(ORG_STRUC.atomType)
                            if ORG_STRUC.STDMOL(i).types(j)==k
                                comp(i,k)=comp(i,k)+1;
                            end
                        end
                    end
                end
                
            end
            %           if ORG_STRUC.varcomp==0
            %                 if ORG_STRUC.molecule==1
            %                     if ~sameComposition(numIons, ORG_STRUC.numIons*comp)
            %                         warningStr = ['Seeds : The number of atoms is inconsistent when reading Seeds-' num2str(loops_seed+1),...
            %                             '\nSkip the rest seeds.'];
            %                         USPEXmessage(0, warningStr, 0);
            %                         return;
            %                     else
            %                         numMols = ORG_STRUC.numIons;
            %                     end
            %                 elseif ~sameComposition(numIons, ORG_STRUC.numIons)
            %                     warningStr =  ['Seeds : The number of atoms is inconsistent when reading Seeds-' num2str(loops_seed+1),...
            %                         '\nSkip the rest seeds.'];
            %                     USPEXmessage(0, warningStr, 0);
            %                     return;
            %                 end
            %             else
            if ORG_STRUC.molecule==1
                if ~isempty(numMolsRead) & length(numMolsRead) ==size( comp, 1 ) & norm(numIons - numMolsRead*comp ) < 1e-3
                    numMols = numMolsRead;
                else
                    numMols = round(numIons/comp);
                end
                numBlocks1 = numMols/ORG_STRUC.numIons;
            else
                numMols = 1;   % fix the bug,  mean nothing here.
                numBlocks1 = numIons/ORG_STRUC.numIons;
            end
            if abs( sum(numBlocks1-round(numBlocks1)) ) > 0.01 | ~isempty(find(numMols<-0.1)) | ~isempty(find(numBlocks1<-0.1))
                warningStr = (['Seeds : Impossible to make combination [ ' num2str(numBlocks1) ' ] when reading Seeds-' num2str(loops_seed+1),...
                    ', so we skip this Seed ...']);
                USPEXmessage(0, warningStr, 0);
                natom = sum(ntyp); % number of atoms
                tmp = fgetl(fid);
                candidate_pop_tmp = fscanf(fid,'%g',[3,natom]);
                tmp = fgetl(fid);  % Read the end the line
                
                %if fgetl(fid)==-1 break; end
                
                continue;
            end
            numBlocks = round(numBlocks1);
            %           end
            
            
            %-------------------------%
            %-- To get atom number  --%
            %-------------------------%            
            natom = sum(ntyp); % number of atoms
            
            %-------------------------%
            %--To get atomic position-%
            %-------------------------%             
            tmp   = fgetl(fid);
            candidate_pop_tmp = fscanf(fid,'%g',[3,natom]);
            candidate_pop_tmp = candidate_pop_tmp';
            tmp  =fgetl(fid); % Read the end the line
            
            %--- maping the atomic positions
            candidate_pop = zeros(size(candidate_pop_tmp));
            startIons=0;
            for i = 1:length(ORG_STRUC.atomType)
                for j = 1:length(atomType_Seeds)
                    if ORG_STRUC.atomType(i)==atomType_Seeds(j)
                        if j==1
                            candidate_pop(startIons+1:startIons+ntyp(j),:)= candidate_pop_tmp(1:ntyp(1),:);
                            startIons = startIons+ntyp(j);
                        else
                            candidate_pop(startIons+1:startIons+ntyp(j),:)= candidate_pop_tmp( sum(ntyp(1:j-1))+1:sum(ntyp(1:j)),: );
                            startIons = startIons+ntyp(j);
                        end
                    end
                end
            end
            
            %--- end of maping the atomic positions
            
            Add = length(POP_STRUC.POPULATION) + 1;
            if POP_STRUC.generation == 1
                ORG_STRUC.initialPopSize = ORG_STRUC.initialPopSize + 1;
                Add = ORG_STRUC.initialPopSize;
            end
            info_parents = struct('parent',{}, 'enthalpy', {});
            info_parents(1).parent = [];
            POP_STRUC.POPULATION(Add).Parents = info_parents;
            POP_STRUC.POPULATION(Add).LATTICE = candidate_lat;
            POP_STRUC.POPULATION(Add).COORDINATES = candidate_pop;
            POP_STRUC.POPULATION(Add).numIons = numIons;
            POP_STRUC.POPULATION(Add).howCome = '  Seeds';
            if ORG_STRUC.molecule==1
                POP_STRUC.POPULATION(Add).numMols = numMols;
                [type, MtypeLIST, numIons]=GetPOP_MOL(numMols);
                POP_STRUC.POPULATION(Add).typesAList = type;
                POP_STRUC.POPULATION(Add).MtypeLIST = MtypeLIST;
                readMOL(Add, 0);
            end
            if ORG_STRUC.varcomp==1
                POP_STRUC.POPULATION(Add).numBlocks = numBlocks;
            end
            loops_seed = loops_seed + 1;
            disp(['seed number ' num2str(loops_seed) ' has been successfully added']);
            
        catch
            warningStr = (['Seeds : Meet a problem when reading Seeds-' num2str(loops_seed+1),...
                ', so we stop the rest seeds reading...']);
            USPEXmessage(0, warningStr, 0);
            break
        end
    end  %while
    status = fclose(fid);
    disp(' ')
    disp('End of pickup Seeds')
    
    [a,b]=unix([ 'echo Generation:' num2str(POP_STRUC.generation) '>> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
    [a,b]=unix([ 'cat POSCARS >> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
    [a,b]=unix([ 'mv POSCARS POSCARS_' num2str(POP_STRUC.generation) ]);
end

cd ..


%
%
%

function numMolsRead = numMolReadSeeds( string0 )


[nothing, numIons] = unix([' echo "' string0 '" |cut -d[  -f2 | cut -d] -f1  ']);
numMolsRead = str2num(numIons);
