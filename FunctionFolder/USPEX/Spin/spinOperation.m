function spinOperation(ini_Ind)


global POP_STRUC
global ORG_STRUC
global OFF_STRUC
global USPEX_STRUC


if ORG_STRUC.spin ~=1
   return;
end

for i = 1:ini_Ind
    if ~isempty(findstr( OFF_STRUC.POPULATION(i).howCome, 'Random') )
       OFF_STRUC.POPULATION(i) = individual_Spin_Init(  OFF_STRUC.POPULATION(i) );
    elseif ~isempty(findstr( OFF_STRUC.POPULATION(i).howCome, 'Heredity') )
       parent = str2num(OFF_STRUC.POPULATION(i).Parents.parent);
%       ID = [];
%       for j = 1:length(POP_STRUC.POPULATION)
%           if abs(POP_STRUC.POPULATION(j).Number-parent(1) ) < 1e-3 | abs(POP_STRUC.POPULATION(j).Number-parent(2) ) < 1e-3
%              ID = [ID, j];
%           end
%           if length(ID) == 2; break; end
%       end
%       ID
%      parent
%       POP_STRUC.POPULATION( ID(1) )
%       POP_STRUC.POPULATION( ID(2) )

       magType1 = USPEX_STRUC.POPULATION( parent(1) ).magmom_ions(end, 1);
       magType2 = USPEX_STRUC.POPULATION( parent(2) ).magmom_ions(end, 1);
       magRatio_Orig = ORG_STRUC.magRatio;
       ORG_STRUC.magRatio = newMagRatio(magType1,magType2);
       OFF_STRUC.POPULATION(i) = individual_Spin_Init(  OFF_STRUC.POPULATION(i) );
       ORG_STRUC.magRatio = magRatio_Orig;
    else
       parent = str2num(OFF_STRUC.POPULATION(i).Parents.parent);
%       ID = [];
%       for j = 1:length(POP_STRUC.POPULATION)
%           if abs(POP_STRUC.POPULATION(j).Number-parent(1) ) < 1e-3
%              ID = j ;
%              break
%           end
%       end
       OFF_STRUC.POPULATION(i).magmom_ions      = USPEX_STRUC.POPULATION( parent ).magmom_ions;
       OFF_STRUC.POPULATION(i).magmom_ions(:,:) = 0;
       OFF_STRUC.POPULATION(i).magmom_ini       = USPEX_STRUC.POPULATION( parent ).magmom_ions(end, :);
       OFF_STRUC.POPULATION(i).magmom_ions(1,:) = OFF_STRUC.POPULATION(i).magmom_ini;
    end
end


Operation = {'Spinmutation'};
Num_Opera = ORG_STRUC.howManySpinmutations;
count = ini_Ind;
for i = 1 : length(Num_Opera)
    for j = 1:Num_Opera(i)
        count = count + 1;
        eval([Operation{i} '(' num2str(count) ')']);
    end
end



%%----------
function magRatio = newMagRatio(magType1,magType2)


magRatio = zeros(1,7);

for magType = [magType1,magType2]

   switch magType
   case  1
      magRatio(1) = magRatio(1)+1;
   case -2
      magRatio(2) = magRatio(2)+1;
   case  2
      magRatio(3) = magRatio(3)+1;
   case -3
      magRatio(4) = magRatio(4)+1;
   case  3
      magRatio(5) = magRatio(5)+1;
   case  4
      magRatio(6) = magRatio(6)+1;
   case  5
      magRatio(7) = magRatio(7)+1;
   end
end


magRatio = magRatio/sum(magRatio);

