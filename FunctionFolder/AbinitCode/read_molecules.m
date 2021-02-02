function read_molecules()
global ORG_STRUC

uff = struct('molecule',{}, 'ZMATRIX', {}, 'name',{}, 'format',{},'radi',{},'types',{}, 'charge', {},...
             'optFlags',{},'index',{}, 'symbol', {}, 'molcenter', {}, 'length', {}, 'CN', {}, ...
             'num_optFlags', {}, 'flex_dihedral', {}, 'bound', {});        
if isempty(ORG_STRUC.numMols)  %molecule generation for 300
   [nothing, results] = unix('ls MOL_[1-9] |wc');
   tmp = str2num(results);
   N_Mols = tmp(1);
   look = 1;
else
   N_Mols = size(ORG_STRUC.numMols,2);
   look = 0; %do not read numMols   
end

for runner =1:N_Mols
    if exist(['MOL_' num2str(runner)])
        handle = fopen(['MOL_' num2str(runner)]);
        name = fgetl(handle);
        if look == 1
           start  = findstr(name, '[');
           ending = findstr(name, ']');
           ORG_STRUC.numMols(runner) = str2num(name(start:ending));
        end
        %how many columns
        if ~isempty(findstr(name, 'charge')) %GULP has one more line to specify charge
           disp('This calculation uses charge model, currently only supported in GULP');
           N_col=9;
        else
           N_col=8;
        end
        if ORG_STRUC.abinitioCode(1)==12 %TINKER
           disp('TINKER needs atom type used in force field, outputs AuxiliaryFiles/#.cssr');
           N_col=9;
        end


        a = fgetl(handle);
        num = sscanf(a,'%*s %*s %*s %g',[1,1]);
        temper=zeros(num,N_col);
        for i=1:num
           a = fgetl(handle);
           if N_col==8
              tmp = sscanf(a,'%s %g %g %g %g %g %g %g');
           elseif N_col==9
              tmp = sscanf(a,'%s %g %g %g %g %g %g %g %g');
           end
           
           type ='';
           shift = -1; %
           if ~isempty(findstr(a,'_'))  %For GULP and DMACRYS
               mark  = findstr(a,'_') - 1; %how many chars for the element
               shift = shift + 2; %_1 has two digits
               uff(runner).symbol(i) = char(tmp(mark+2));
           else 
               mark = length(tmp) + 1 - N_col;           
           end

           for ind=1:mark
               type(ind)=char(tmp(ind));
               shift=shift+1;
           end

           temper(i,1) = elementID(type);  %%to find the index for the element

           for k= 2:N_col
               temper(i,k)=tmp(k+shift);
           end

           if ORG_STRUC.abinitioCode(1)==6  %DMACRYS
              uff(runner).index(i) = i;
           end

        end

        uff(runner).name = name;
        types = zeros(num,1);
        for ind = 1:num
            types(ind)= find(ORG_STRUC.atomType==temper(ind,1));
        end
        uff(runner).types = types;
        uff(runner).molecule = temper(:,2:4);
        uff(runner).format = temper(:,5:7);
        uff(runner).molcenter = mean(uff(runner).molecule);
    
        if N_col==9  %GULP+CHARGE
           uff(runner).charge = temper(:,9); %Charge for GULP, type for TINKER
        else
           uff(runner).charge = zeros(i,1);
        end


        optFlags = zeros(num,3);
        if ORG_STRUC.dimension==1
           uff(runner).bound = temper(:,8);
        else   %if ORG_STRUC.abinitioCode(1)~=6
           optFlags(:,3) = temper(:,8);
        end
        optFlags(1,:) = [1 1 1];
        if num==2
           optFlags(2,:) = [0 1 1];
        elseif num>=3
           optFlags(3,:) = [0 0 1];
        end
        uff(runner).optFlags = optFlags;
        uff(runner).num_optFlags = sum(sum(optFlags)) - 6;
        mark = 0;
        if uff(runner).num_optFlags > 0
           for i=4:num
               for j=1:3
                   if uff(runner).optFlags(i,j)==1
                      mark = mark +1 ;
                      uff(runner).flex_dihedral(mark,1:2)=[i,j];
                   end
               end
           end
        end


        %%%%%%%%%%%%%% HERE THE DISTANCES NEED TO BE CALCULATED

       distances = zeros(num,1);
       for ind = 1:num
           distances(ind)=norm(temper(ind,2:4)-uff(runner).molcenter);
       end
       if num==1
           uff(runner).radi = 0.8*str2num(covalentRadius(ORG_STRUC.atomType(types)));
       else
           uff(runner).radi = max(distances)+0.4;
       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%   Zmatrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (ORG_STRUC.dimension==1) & (runner==1)
          uff(runner) = ReviseMOL(uff(runner), ORG_STRUC.atomType);
       elseif num > 1 %orientate the molecule according to the priciple Rotation axis
          MOLCOORS = uff(runner).molecule;
          MOLCCORS = bsxfun(@minus, MOLCOORS, mean(MOLCOORS)); % Matrix-
          [a,b] = PrincipleAxis(MOLCOORS);
          uff(runner).molecule = MOLCOORS*a;
       end

       format = uff(runner).format;
       coords = uff(runner).molecule;
       types  = uff(runner).types;
       uff(runner).ZMATRIX = real(NEW_coord2Zmatrix(coords, format));
       disp('Please see the calculated Zmatrix in OUTPUT.txt')
       radiu = zeros(num,1);
       for i = 1:num
           radiu(i) = str2num(covalentRadius(ORG_STRUC.atomType(types(i))));
       end
       uff(runner).CN = find_pair(coords, radiu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       status = fclose(handle);

    else
        disp(['For Molecules, you did not give enough MOL_x files']);
        disp(['The program is stopping.........']);
        [nothing, nothing] = unix(['echo For Molecules, you did not give enough MOL_x files >>' ORG_STRUC.resFolder '/' ORG_STRUC.log_file]);
        [nothing, nothing] = unix(['echo The program is stopping......... >>' ORG_STRUC.resFolder '/' ORG_STRUC.log_file]);
        quit
    end
end

ORG_STRUC.STDMOL=uff;

for i=1:N_Mols
    for j=i:N_Mols
        ORG_STRUC.CenterminDistMatrice(i,j) = uff(i).radi + uff(j).radi;
        ORG_STRUC.CenterminDistMatrice(j,i) = uff(i).radi + uff(j).radi;        
    end
end
try
    %[nothing, CenterminDistMatrice] = unix(['./getFromTo MolCenters EndMol INPUT.txt']);
    CenterminDistMatrice=python_uspex('FunctionFolder/getInput.py', ['-f INPUT.txt -b MolCenters -e EndMol']);
    if ~isempty(CenterminDistMatrice)
      if size(ORG_STRUC.CenterminDistMatrice) ~= [N_Mols, N_Mols]
        disp(['MolCenters tag has wrong dimension, Program STOPS']);
        quit
      else
        ORG_STRUC.CenterminDistMatrice = str2num(CenterminDistMatrice);
      end
    end
catch
end

