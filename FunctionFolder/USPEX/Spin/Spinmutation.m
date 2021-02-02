function spinMutation(Ind_No)

% USPEX Version 7.4.1

global POP_STRUC
global ORG_STRUC
global OFF_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING Mutants by atom positions mutation%%%%%%%%%%%%%%%%%


goodMutSpin = 0;
safeguard = 0;


% magType : 
%      NM:  1, 
%   FM-LS: -2, 
%   FM-HS:  2, 
%  AFM-LS: -3, 
%  AMF-HS:  3, 
% FM-HSLS:  4, 
% AFM-HSLS: 5,

mutTryStep = 1;
while goodMutSpin ~= 1

      toMutate = find(ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
      ind     = POP_STRUC.ranking(toMutate(end));
      numIons = POP_STRUC.POPULATION(ind).numIons;
      magmom_Mut = zeros(1,1+sum(numIons));
      magmom     = POP_STRUC.POPULATION(ind).magmom_ions(end,:);

      magTypeTrue=magTypeIdentify(magmom, POP_STRUC.POPULATION(ind).magmom_ini);

      mutTryStep=0;
      
      while ~goodMutSpin 
           magRatio    = ORG_STRUC.magRatio;
           switch magTypeTrue(1)
           case  1
              magRatio(1) = 0;
           case -2
              magRatio(2) = 0;
           case  2
              magRatio(3) = 0;
           %case  3 or -3
               % AFM can mutate to AFM
           case  4  
              magRatio(6) = 0; 
           case  5
              magRatio(7) = 0; 
           end 
           magmom_Mut  = initialize_magMom( numIons, magRatio );
           spinMutType=[magTypeString(magTypeTrue(1)) '->' magTypeString(magmom_Mut(1))];

           if abs(magmom_Mut(1))==3 & abs(magTypeTrue(1))==3
                goodMutSpin = checkAFMmutat(magmom_Mut, magTypeTrue);         
           else
                goodMutSpin = 1;
           end 

           if mutTryStep>10 
              % to avoid very small ratio of magnetic type, for example:
              % [ 0 1 1 0 0 0.01 0]
              disp('Warning: exceed to the max spin muataion try!')
              goodMutSpin=1;
              spinMutType=[magTypeString(magTypeTrue(1)) '->' magTypeString(magmom_Mut(1))];
           end
           mutTryStep = mutTryStep + 1;
      end 

    if goodMutSpin == 1
        OFF_STRUC.POPULATION(Ind_No).COORDINATES = POP_STRUC.POPULATION(ind).COORDINATES;
        OFF_STRUC.POPULATION(Ind_No).LATTICE     = POP_STRUC.POPULATION(ind).LATTICE; 
        OFF_STRUC.POPULATION(Ind_No).numIons     = POP_STRUC.POPULATION(ind).numIons;

        try
          OFF_STRUC.POPULATION(Ind_No).numBlocks   =POP_STRUC.POPULATION(ind).numBlocks;
        end


        OFF_STRUC.POPULATION(Ind_No).howCome = 'Spinmutate';
        OFF_STRUC.POPULATION(Ind_No).magmom_ini = magmom_Mut;
        OFF_STRUC.POPULATION(Ind_No).magmom_ions(1,:) = magmom_Mut;

        info_parents = struct('parent', {},'spinMutType', {}, 'enthalpy', {});
        info_parents(1).parent   = num2str(POP_STRUC.POPULATION(ind).Number);
        info_parents.spinMutType = spinMutType;
        info_parents.enthalpy = POP_STRUC.POPULATION(ind).Enthalpies(end)/sum(numIons);
        OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
        disp(['Structure ' num2str(Ind_No) ' generated from Spinmutation : ' spinMutType ]);

    end

end

%%%%%%%%%%%%%%% END creating mutants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isOK = checkAFMmutat(orgMag, mutMag)

seq=find( abs(orgMag(2:end))>0 );

if orgMag(2)<0   %% To make sure the 1st atom have a spin up magmom
   orgMag(2:end)=-orgMag(2:end);
end

orgMag(seq+1)=orgMag(1+seq)./abs( orgMag(1+seq) );
mutMag(seq+1)=mutMag(1+seq)./abs( mutMag(1+seq) );


if sum(abs(orgMag(2:end)-mutMag(2:end)))~=0
   isOK = 1;
else
   isOK = 0;
end



