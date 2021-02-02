function checkCompatibility()


%
% This function is used to keep the compatibility of different USPEX versions.
%
% created by Guangrui Qian
%

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC



% svn. 397
if ~isfield( ORG_STRUC, 'lattice')
    ORG_STRUC.lattice = ORG_STRUC.constLat;
end

% svn. 557
load( [POP_STRUC.resFolder '/USPEX.mat'] );
if ~isfield( USPEX_STRUC, 'SYSTEM')
	atomType = ORG_STRUC.atomType;
	USPEX_STRUC.SYSTEM(1).atomType = atomType;
	for i = 1:length(atomType)
	   USPEX_STRUC.SYSTEM(1).atomType_symbol{i} = megaDoof(atomType(i));
	end
	USPEX_STRUC.SYSTEM(1).Fp_weight = ORG_STRUC.weight;
 	safesave([POP_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
end

% fixRndSeed feathure compatibility
if ~isfield( ORG_STRUC, 'fixRndSeed')
	ORG_STRUC.fixRndSeed = 0;
end

% firstGeneSplitAll
if ~isfield( ORG_STRUC, 'firstGeneSplitAll')
    N_T = size(ORG_STRUC.numIons,1);
    splitting = zeros(1,N_T);
    findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
end

% COPEX plugin compatibility
if ~isfield(ORG_STRUC, 'pluginType')
    ORG_STRUC.pluginType = 0;
end
if ~isfield(ORG_STRUC, 'currentGenDone')
    ORG_STRUC.currentGenDone = 0;
end
if ~isfield(ORG_STRUC, 'startNextGen')
    ORG_STRUC.startNextGen = 1;
end

% svn.664 protein compatibility
if ~isfield(ORG_STRUC, 'howManySecSwitch')
    ORG_STRUC.howManySecSwitch = 0;
end
if ~isfield(ORG_STRUC, 'howManyShiftBorder')
    ORG_STRUC.howManyShiftBorder = 0;
end
if ~isfield(ORG_STRUC, 'fracSecSwitch')
    ORG_STRUC.fracSecSwitch = 0;
end
if ~isfield(ORG_STRUC, 'fracShiftBorder')
    ORG_STRUC.fracShiftBorder = 0;
end

% spin  feathure compatibility
if ~isfield( ORG_STRUC, 'spin')
    ORG_STRUC.spin = 0;
end  
if ~isfield( ORG_STRUC, 'fracSpin')
    ORG_STRUC.fracSpin = 0;
end
if ~isfield( ORG_STRUC, 'howManySpinmutations')
    ORG_STRUC.howManySpinmutations = 0;
end
if ~isfield( ORG_STRUC, 'ldaU')
    ORG_STRUC.ldaU = 0;
end 

% collct forces compatibility 
if ~isfield( ORG_STRUC, 'collectForces')
   ORG_STRUC.collectForces = 0;
end


% magnetci feature compatibility
%if  ~isfield( POP_STRUC.POPULATION(1), 'magmon_ini')
%    for i = 1:length(POP_STRUC.POPULATION)
%        POP_STRUC.POPULATION(i).magmon_ini = [];
%        POP_STRUC.POPULATION(i).magmon_mom = [];
%    end
%end

%%
%% Do we need to save ORG_STRUC and POP_STRUC to the .mat file here?
%%
