function magmomTrue = magTypeIdentify(magmom, mag_ini)


%
%  This fucntion is used to identify what type magnetic system is
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  magmom_ions format:[magType, X X X X 0 0 0 0 ....]
%
% magType : 
%      NM:  1, 
%   FM-LS: -2, 
%   FM-HS:  2, 
%  AFM-LS: -3, 
%  AMF-HS:  3, 
% FM-HSLS:  4, 
% AFM-HSLS: 5,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


magTypeOrig = magmom(1);
magTypeTrue = 0;
magmomTrue  = zeros(1,length(magmom));

criteriaLSHS=1.5;


seqMagAtom     = ( find(abs(mag_ini(2:end))>0) + 1);
if isempty(seqMagAtom) 
   seqMagAtom  = ( find(abs(magmom(2:end))>0) + 1);
end
numMagAtom=length( seqMagAtom );

if 1==seqMagAtom
       magmomTrue(1)=1;
       magmomTrue(2:end)=0;   %NM   
   return;
end

%%----------------------
meanMagMom = mean(abs(magmom(seqMagAtom)));
meanMagErr = mean( abs( abs(magmom(seqMagAtom))-meanMagMom ) );
isnotAFM   = sum(magmom(seqMagAtom)./abs(magmom(seqMagAtom)));

if  ( sum(abs(magmom(2:end)))<0.03*numMagAtom ) | isempty(seqMagAtom)
    magTypeTrue=1;
    magmomTrue(1)=1;
    magmomTrue(2:end)=0;   %NM
elseif  (abs( sum(magmom(2:end)) ) < 0.25*numMagAtom) & (abs(isnotAFM)<1e-2)
    if      max(abs(magmom(seqMagAtom)))<criteriaLSHS
           magTypeTrue=-3;     %AFM-LS
    elseif  min(abs(magmom(seqMagAtom)))>criteriaLSHS
           magTypeTrue= 3;     %AFM-HS
    else
           magTypeTrue= 5;     %AFM-HSLS
    end
  else 
   %magmom(seqMagAtom)
   %mag_ini
    if     max(abs(magmom(seqMagAtom)))<criteriaLSHS
           magTypeTrue=-2;     %FM-LS
    elseif min(abs(magmom(seqMagAtom)))>criteriaLSHS
           magTypeTrue= 2;      %FM-HS
    else
           magTypeTrue= 4;     %FM-HSLS
    end
end

if magmom(2) < 0
    magmom(2:end) = -magmom(2:end);
end
magmomTrue = [ magTypeTrue, magmom(1,2:end) ];

magString=magTypeString(magTypeTrue);
%disp([magString,' : ', num2str(magmomTrue)])
