function loops_seed = Pick_Seeds_M200()

% USPEX Version 9.3
% Change: block varcomp support added
% Change: the seeds structure would be added into the end of the POPULATION
% Change: VASP-5.2 format 
global POP_STRUC
global ORG_STRUC
warning off

cd Seeds
loops_seed = 0;

if exist('POSCARS')
   [fid,message] = fopen('POSCARS');
   while 1
     try
         tmp = fgetl(fid); % system description
         scale_factor = fgetl(fid); % 1 = numbers in angstrems
         optlattice = fscanf(fid,'%g',[3,3]); % lattice vectors
         candidate_lat = optlattice'*str2num(scale_factor);
         tmp = fgetl(fid);  % we don't need this line
         atomType = fgetl(fid); % VASP5, elemental type
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
             warningStr = ['Seeds : Elemental types are inconsistent in Seeds-' num2str(loops_seed+1) ];
             USPEXmessage(0, warningStr, 0);
             break;
         end

%%%%%    End AtomType     
         numIons = zeros(1,length(ORG_STRUC.atomType));
         for i=1:length(ORG_STRUC.atomType)
             for j=1:length(ntyp)
                 if atomType1(j)==ORG_STRUC.atomType(i)
                     numIons(i) = ntyp(j);
                 end 
             end
         end

         if numIons ~= ORG_STRUC.numIons
            warningStr = ['NumIons are inconsistent in Seeds-' num2str(loops_seed+1) ];
            USPEXmessage(0, warningStr, 0);
            break;
         end
%%%%%    End numIons    

         natom = sum(ntyp); % number of atoms
         tmp = fgetl(fid);
         candidate_pop = fscanf(fid,'%g',[3,natom]);
         candidate_pop = candidate_pop';
         tmp = fgetl(fid);
         Add = length(POP_STRUC.POPULATION) + 1;
         if POP_STRUC.generation == 1
            ORG_STRUC.initialPopSize = ORG_STRUC.initialPopSize + 1;
         end
         info_parents = struct('Seeds',{}, 'enthalpy', {});
         loops_seed = loops_seed + 1;
         info_parents(1).Seeds = loops_seed;
         POP_STRUC.POPULATION(Add).Parents = info_parents;
         POP_STRUC.POPULATION(Add).LATTICE = candidate_lat;
         POP_STRUC.POPULATION(Add).COORDINATES = candidate_pop;
         POP_STRUC.POPULATION(Add).numIons = numIons;
         disp(['seed number ' num2str(loops_seed) ' has been successfully added']);

     catch
        break
     end
   end  %while
   fclose(fid); 
   disp('end of pickup Seeds')
   [nothing, nothing] = unix([ 'echo Generation:' num2str(POP_STRUC.generation) '>> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
   [nothing, nothing] = unix([ 'cat POSCARS >> ../' ORG_STRUC.resFolder '/Seeds_history' ]);
   [nothing, nothing] = unix([ 'rm POSCARS' ]);
end

cd ..
