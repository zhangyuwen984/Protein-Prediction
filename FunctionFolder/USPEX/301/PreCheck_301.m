function PreCheck_301()

% This function is used to check the INPUT.txt parameters and some specific
% conditions 


global ORG_STRUC


%-- We no longer support single block in 301

if size(ORG_STRUC.numIons,1) == 1
    error('==> USPEX no longer supports single block in 301, please use 300 instead! <==');
end

