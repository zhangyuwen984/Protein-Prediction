function ShiftBorder(Ind_No)
%function ShiftBorder()
%Ind_No = 6
%Ind_No

% $Rev: 673 $
% $Author: maxim $
% $Date: 2014-10-22 02:42:57 +0400 (Wed, 22 Oct 2014) $

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CREATING Offspring with Secondary Switch %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angles_num = size(POP_STRUC.POPULATION(1).ANGLES(:,:), 1);

% Select a structure using tournament scheme:
toMutate = find(ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
ind      = POP_STRUC.ranking(toMutate(end));

sec_struct     = POP_STRUC.POPULATION(ind).SEC_STRUCT;
rand_structure = POP_STRUC.POPULATION(ind).ANGLES(1:end,:);

% Find continuous secondary structures: 
cont_seq = struct('names'  , '', ...
                  'numbers', [] ...
                  );
start_s = 1;
counter = 1;  %Cycle to determine secondary structures boundaries in the chain
for i=1:size(sec_struct, 2)
    if i > 1
        if ~isequal(sec_struct{i}, sec_struct{i-1}) %if new structure was found, we write down in the cont_seq struc length and name of the previous
            end_s = i-1;
            diap = [start_s, end_s];
            start_s = i;
            cont_seq.names{counter} = sec_struct{i-1};
            cont_seq.numbers = [cont_seq.numbers; diap];
            counter = counter + 1;
        end
        if i == size(sec_struct, 2) %the special variant for the last aminoacid
            end_s = i;
            diap = [start_s, end_s];
            cont_seq.names{counter} = sec_struct{i};
            cont_seq.numbers = [cont_seq.numbers; diap];
            counter = counter + 1;
        end
    end
end

potentialOffspring = rand_structure;

if size(cont_seq.numbers) ~= [0, 0]

    %disp(['Secondary structures (each residue): ' strjoin(sec_struct, ' ')]);
    %disp(['Secondary structures (compressed)  : ' strjoin(cont_seq.names, ' ')]);

    disp(['Secondary structures (each residue): ' strrep(reshape(char(sec_struct)',1,[]), ' ',' ')    ]);
    disp(['Secondary structures (compressed)  : ' strrep(reshape(char(cont_seq.names)',1,[]), ' ',' ')]);

    if size(cont_seq.names,2) > 1
        borders_num       = size(cont_seq.names,2) - 1; %number of boundaries
        borders2shift_max = round(borders_num/2); %number of boundaries to shift
        borders2shift_num = RandInt(1,1,[1,borders2shift_max]); %randomly chose the number of shifted boundaries
        borders2shift     = randperm(borders_num, borders2shift_num); %which borders will be shift
        plus_minus        = RandInt(1,size(borders2shift,2),[1,2]); %randomly chose the direction of the shift and number of residues
        plus_minus        = (plus_minus - 1.5) * 10

        disp(['-- We have ' num2str(borders_num) ' borders. Total ' num2str(borders2shift_num) ' border(s) will be shifted:']);
        for i=1:borders2shift_num
            border_id = borders2shift(i);
            left_right = plus_minus(i);
            if left_right < 0
                move_direction = 'left';
            elseif left_right > 0
                move_direction = 'right';
            end
            if left_right < 0     % move border to the left
                sec2change = border_id + 1;
                move_id    = cont_seq.numbers(sec2change,1);    % first pair of angles of the right secondary structure
                next = -1;
                if abs(cont_seq.numbers(sec2change-1,1)-cont_seq.numbers(sec2change,1)) < abs(plus_minus(i))
                    left_right = cont_seq.numbers(sec2change-1,1)-cont_seq.numbers(sec2change,1);
                end
            elseif left_right > 0  % move border to the right
                sec2change = border_id;
                move_id    = cont_seq.numbers(sec2change,2);    % last pair of angles of the left secondary structure
                next = 1;
                if abs(cont_seq.numbers(sec2change+1,2)-cont_seq.numbers(sec2change,2)) < abs(plus_minus(i))
                    left_right = cont_seq.numbers(sec2change+1,2)-cont_seq.numbers(sec2change,2);
                end
            end
            disp(['---- Border #' num2str(border_id) ' will be moved on ' num2str(abs(left_right)) ' residues to the ' move_direction ': ' cont_seq.names{sec2change} ' -> ' cont_seq.names{sec2change+next}]);
            for k=1:abs(left_right)
                 potentialOffspring(move_id + sign(left_right) * k,:) = potentialOffspring(move_id,:);       
            end
        end
        
    else
        disp('Shift Border cannot be applied since only one secondary structure was found. Execute SecSwitch operator.')
        SecSwitch(Ind_No);
        return;
    end
else
    disp('Shift Border was not applied: no secondary structures found by STRIDE.');
end

OFF_STRUC.POPULATION(Ind_No).ANGLES   = potentialOffspring;
% We need just first list of amino acids, since it won't change:
OFF_STRUC.POPULATION(Ind_No).RESIDUES = POP_STRUC.POPULATION(1).RESIDUES;
OFF_STRUC.POPULATION(Ind_No).numIons  = ORG_STRUC.numIons;
OFF_STRUC.POPULATION(Ind_No).howCome  = 'ShiftBorder';

info_parents = struct('parent', {}, 'enthalpy', {});
info_parents(1).parent = num2str(POP_STRUC.POPULATION(ind).Number);
info_parents.enthalpy = POP_STRUC.POPULATION(ind).Enthalpies(end)/ORG_STRUC.numIons;
OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;


disp(['Structure ' num2str(Ind_No) ' generated by Shift Border']);
