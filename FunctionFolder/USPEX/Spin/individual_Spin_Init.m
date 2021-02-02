function POPULATION = individual_Spin_Init( POPULATION ) 



global ORG_STRUC

%%%%%%%%%%%%%%%%%%%% START initial SPIN MagMON %%%%%%%%%%%%%%%%%
% RANDSEED for MAGMOM TYPE !! Condition:
%
% 0  < randSeed <= 0.2 NM :  MAGMOM= 0...
% 0.2< randSeed <= 0.4 FM-LS :  MGGMOM= 1 1 1 ...  
% 0.4< randSeed <= 0.6 FM-HS :  MGGMOM= 4 4 4 ...
% 0.6< randSeed <= 0.8 AFM-LS :  MGGMOM= -1 1 -1 1 ...
% 0.8< randSeed <= 1   AFM-HS :  MGGMOM= -4 4 -4 4 ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  magmom_ions format
%
%     N(ions)   s  p  d  tot
%%---------------------------------------
%disp('Here is the magnetic moments initialization!!')


POPULATION.magmom_ini  = zeros( 1, 1+sum(POPULATION.numIons) );  
POPULATION.magmom_ions = zeros( length([ORG_STRUC.abinitioCode]), 1+sum(POPULATION.numIons) );  
POPULATION.magmom_ions(1,:) = initialize_magMom( POPULATION.numIons, ORG_STRUC.magRatio ); 
POPULATION.magmom_ini       = POPULATION.magmom_ions(1,:);

POPULATION.mag_moment       = zeros( 1,length([ORG_STRUC.abinitioCode]) );

%%%%%%%%%%%%%%%%%%%% END SPIN MagMON %%%%%%%%%%%%%%%%%

