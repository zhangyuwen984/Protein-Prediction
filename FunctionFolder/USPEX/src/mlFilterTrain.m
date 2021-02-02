function mlFilterTrain()

global ORG_STRUC;
global POP_STRUC;
    
assert(ORG_STRUC.mlPeerFilter > 0, ...
       'mlFilterTrain should NOT be called as mlPeerFilter is not set.');

assert(ORG_STRUC.mlPeerFilter < 2, ...
       'mlPeerFilter should be set to 1.');

model_names = {'krr'};  % TODO: impl. models 'knr'; 'svr'; 'nn'
model = model_names{ORG_STRUC.mlPeerFilter};
    
model_file = [ORG_STRUC.resFolder '/uspexml-' model '.pkl'];
uspex_mat = [ORG_STRUC.resFolder '/USPEX.mat'];
disp(['Training a ' model ' model with structures before generation ' ...
		    num2str(POP_STRUC.generation) '.']);
train_cmd = ['uspexml train -m ' model ' -f ' model_file ...
				 ' -d ' uspex_mat ' -v 5 -i 128 -l 8'];
[s, w] = unix(train_cmd);
if 0 ~= s
  disp(['Error running ' train_cmd]);
end
    
