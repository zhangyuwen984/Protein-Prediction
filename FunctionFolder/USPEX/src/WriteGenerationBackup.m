function WriteGenerationBackup()

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC
global POOL_STRUC

safesave ([POP_STRUC.current_dir '/POP_STRUC.mat'],  POP_STRUC)
safesave ([POP_STRUC.current_dir '/ORG_STRUC.mat'],  ORG_STRUC)
safesave ([POP_STRUC.current_dir '/USPEX.mat'],    USPEX_STRUC)
safesave ([POP_STRUC.current_dir '/POOL.mat'],      POOL_STRUC)
