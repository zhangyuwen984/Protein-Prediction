function magmom_ions = initialize_Magmom( numIons, magRatio_ORIG)

% USPEX Version 9.3.6
% 

global ORG_STRUC

%---------------------------------------------- 
% RANDSEED for MAGMOM TYPE !! Condition:
%
% 0  < randSeed <= 0.1     NM :  MAGMOM=  0 ... 
% 0.1< randSeed <= 0.3  FM-LS :  MAGMOM=  1  1   1  1 ...  
% 0.3< randSeed <= 0.5  FM-HS :  MAGMOM=  4  4   4  4 ...
% 0.5< randSeed <= 0.75 AFM-LS :  MAGMOM= -1  1  -1  1 ...
% 0.75< randSeed <= 1   AFM-HS :  MAGMOM= -4  4  -4  4 ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  magmom_ions format:[magType, X X X X 0 0 0 0 ....]
%magType
% NM: 1, FM-LS: -2, FM-HS: 2, AFM-LS: -3, AMF-HS: 3
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% atomType : 21 - 30, 39 - 48 will set Magmom only
% only "ONE" magnetic atoms can be set !!!!

sumIons    = sum(numIons);
atomType   = ORG_STRUC.atomType;

if sum(magRatio_ORIG)< 1e-3
   magRatio_ORIG = ORG_STRUC.magRatio;
end
magRatio_ORIG = magRatio_ORIG/sum(magRatio_ORIG);


randSeq=zeros(1,sumIons); randSeq(:)=1;

nStart=0;
numMagAtom=0;

for i = 1:length(numIons)
    if ( (atomType(i)<21) | (atomType(i)>49) | ((atomType(i)>30) & (atomType(i)<39)) )  &  length(atomType) > 1
        if 1==i
           nStart=1;
        else
           nStart=sum(numIons(1:i-1))+1;
        end
        randSeq(nStart:nStart+numIons(i)-1)=0;
    elseif numIons(i)>0
        numMagAtom=numMagAtom+1;
        if 1==i
           nStart=1;
        else
           nStart=sum(numIons(1:i-1))+1;
        end
        tmp=randperm(numIons(i));
        seq=tmp( floor(numIons(i)/2)+1:numIons(i) )+nStart-1;
        randSeq( seq ) =  -1;
        if randSeq(nStart)== -1
            randSeq(nStart:nStart+numIons(i)-1)=randSeq(nStart:nStart+numIons(i)-1)*(-1);
        end
    end
end

if numMagAtom>1
   disp('More than one atoms have to be set with magnetic, this is not implemented yet!!');
   exit();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
magmom_ions=zeros(1,1+sumIons);

for it = 1:length(magRatio_ORIG)
    magRatio(it)=sum(magRatio_ORIG(1:it));
end
randSeed=rand();
if       (randSeed<=magRatio(1))                    %%NM : non-magm
       magmom_ions(1)=1;
       magmom_ions(2:end)= 0;
elseif   (randSeed<=magRatio(2)) & (randSeed>magRatio(1))   %%FM-LS
       magmom_ions(1)=-2;
       magmom_ions(2:end)= abs(randSeq)*1 ;           
elseif   (randSeed<=magRatio(3)) & (randSeed>magRatio(2))   %%FM-HS
       magmom_ions(1)=2;
       magmom_ions(2:end)= abs(randSeq)*4 ;
elseif   (randSeed<=magRatio(4)) & (randSeed>magRatio(3))   %%AFM-LS
       magmom_ions(1)=-3;
       magmom_ions(2:end)= randSeq*1 ;
       if sum(magmom_ions(2:end)) ~= 0 
          disp(['Warning: odd magnetic atoms cannot be set as AFM-LS type, set as FM-LS.'])
                 magmom_ions(1)=-2;
                 magmom_ions(2:end)= abs( magmom_ions(2:end)  ) ;
       end
elseif   (randSeed<=magRatio(5)) & (randSeed>magRatio(4))   %%AFM-HS
       magmom_ions(1)=3;
       magmom_ions(2:end)= randSeq*4 ;
       if sum(magmom_ions(2:end)) ~= 0 
          disp(['Warning: odd magnetic atoms cannot be set as AFM-HS type, set as FM-HS.'])
                 magmom_ions(1)=2;
                 magmom_ions(2:end)= abs( magmom_ions(2:end)  ) ;
       end
elseif   (randSeed<=magRatio(6)) & (randSeed>magRatio(5))   %%FM-LSHS
       magmom_ions(1)=4;
       magmom_ions(2:end)=randSeq;
       seq= find( abs(randSeq)>0 );
       temp=randperm( sum(abs(randSeq)) ); 
       mutSeq=temp( 1:abs(round(normrnd(0.5,0.1,1,1)*length(seq))) );
       magmom_ions( seq(mutSeq)+1 )= 4 ;
       magmom_ions(2:end)=abs( magmom_ions(2:end) );
elseif   (randSeed<=magRatio(7)) & (randSeed>magRatio(6))   %%AFM-LSHS
       magmom_ions(1)=5;
       magmom_ions(2:end)= randSeq*4 ;
       if sum(magmom_ions(2:end)) ~= 0
          disp(['Warning: odd magnetic atoms cannot be set as AFM-HSLS type, set as FM-HSLS.'])
                 magmom_ions(1)=4;
                 magmom_ions(2:end)= abs( magmom_ions(2:end)  ) ;
       end
end


if sum(abs(magmom_ions(2:end)))==0
    magmom_ions(1)=1;
end

%fprintf(1, [' ', magTypeString(magmom_ions(1)), ' : ', num2str(magmom_ions(2:end)), '\n']);

%fprintf(1, [' ', magStr, '\n']);
