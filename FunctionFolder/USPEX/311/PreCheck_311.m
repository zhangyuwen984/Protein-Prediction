function PreCheck_311()

% This function is used to check the INPUT.txt parameters and some specific
% conditions 


global ORG_STRUC


%-- We no longer support single block in 311

if size(ORG_STRUC.numIons,1) == 1
    error(' ==> USPEX no longer supports single block in 311, please use 310 instead! <==');
end