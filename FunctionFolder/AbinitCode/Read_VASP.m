function target = Read_VASP(flag, ID)
%-100: if VASP is complete (Seems not used any more)
% 0: if 1st SCF is done & SCF is converged
% 1: energy/enthalpy
%-1: collect all energy
% 2: pressue tensor 
%-2: all pressure tensor
% 3: dielectric constant
% 4: gap
% 5: magmoment
%-5: atomic magmoments
% 6: elastic constant matrix (array data)
% 7: atomic forces
% 8: relaxation atomic positions
% 9: relaxation lattice paramter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if     flag == -100
    [nothing, results] = unix('grep Voluntary OUTCAR');
    if isempty(results)
        %          disp('VASP OUTCAR is not completely Done');
        USPEXmessage(1101,'',0);
        [nothing, nothing] = unix(['cp OUTCAR ERROR-OUTCAR-' ID]);
        target = 0;
    else
        target = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag == 0
    [nothing, results1] = unix('./getStuff OSZICAR " F= " 4');
    [nothing, results2] = unix('grep F OSZICAR');
    %Error mess: Binary file OSZICAR matches
    if isempty(results1) | findstr(results2,'Binary file')
        %disp('VASP 1st SCF is not done or OSZICAR/OUTCAR is broken');
        USPEXmessage(1102,'',0);
        [nothing, nothing] = unix(['cp OUTCAR ERROR-OUTCAR-' ID]);
        target = 0;
    else
        [nothing, whichLine] = unix(['grep -n "F=" OSZICAR | tail -n 1  | cut -d":" -f1 ']);
        [nothing, nothing] = unix( ['sed -n " ' num2str((str2num(whichLine)-1)) ' p" OSZICAR > VASPstepTmp']);
        [nothing, vaspSCFsteps] = unix('./getStuff VASPstepTmp ": " 2 | tail -n 1');
        vaspSCFsteps = str2num(vaspSCFsteps);
        [nothing, nothing] = unix('rm VASPstepTmp');
        [nothing, NELM]=unix('./getStuff OUTCAR " NELM " 4');
        NELM = str2num(strrep(NELM, ';', ''));
        if vaspSCFsteps<NELM
            target = 1;
        else
            USPEXmessage(1103,'',0);
            %disp('VASP SCF is not converged');
            [nothing, nothing] = unix(['cp OUTCAR ERROR-OUTCAR-' ID]);
            target = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  1 | flag ==  -1%
    [nothing, Energy_Str] = unix('./getStuff OSZICAR " F= " 4');
    Energy = str2num(Energy_Str);
    if flag == 1
        target = Energy(end);
    else
        target = Energy;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  2 | flag ==  -2% Pressure Tensor
    %[nothing, prestenStr] = unix('awk -f VASP_pres.awk OUTCAR');
    %target = str2num(prestenStr);
    stress = callAWK('VASP_pres.awk','OUTCAR');
    if flag == 2
        target = stress(end-2:end, 1:3);
    else
        target = stress;
    end

    %USPEXmessage(1104,'',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  3 % Dielectric Constant
    target = Read_VASP_Diel();
    %USPEXmessage(1105,'',0);
elseif flag ==  4 % Gap
    gap = Read_VASP_EIGENVAL();
    %USPEXmessage(1106,'',0);
    if gap < 0.1
       [nothing, numStr] = unix('sed -n "6p" POSCAR');
       NCoords= sum( str2num(numStr) );
       [gap_tmp, fermiDOS] = Read_VASP_DOSCAR(NCoords);
       gap = gap - fermiDOS;
    end
    target = gap;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag ==  5 % Mag_moment
    [nothing, magmomStr] = unix('./getStuff OSZICAR mag 11');
    magmom = str2num(magmomStr);
    if isempty(magmom)
        USPEXmessage(1107,'',0);
        target = 0;
    else
        target=magmom(end); 
    end
elseif flag == -5 % Atomic Mag_moment

    % reads MAGMOM in  OUTCAR file. Format:
    % magnetization (x)
    %
    %# of ion     s       p       d       tot
    %----------------------------------------
    %  1        0.008   0.009   2.911   2.927
    %  2        0.008   0.009   2.911   2.927
    %  3       -0.003  -0.003   0.678   0.672
    %  4       -0.003  -0.003   0.678   0.672
    %  5        0.018   0.114   0.000   0.132
    %  6        0.019   0.113   0.000   0.132
    %  7        0.018   0.114   0.000   0.132
    %  8        0.019   0.113   0.000   0.132
    %------------------------------------------------
    %tot        0.082   0.467   7.178   7.727

    [nothing, numStr] = unix('sed -n "6p" POSCAR');
    numIons= sum( str2num(numStr) );

    magStr=callAWK('VASP_mag.awk', 'OUTCAR', ['num=', num2str(numIons)]);
    if isempty(magStr)
        USPEXmessage(1107,'',0);
        target=zeros(numIons,5);
        target(1:end,1)=1:numIons;
    else
        target=magStr(end-numIons+1:end,:);
    end

elseif flag ==  6 % Elastic properties
    elas_orig = callAWK('VASP_elas.awk', 'OUTCAR');   % [nothing, nothing] = unix('awk -f VASP_elas.awk OUTCAR');
    try
        elasticMatrix = elas_orig([1,2,3,5,6,4],:)/10.0;
    catch
        USPEXmessage(1108,'',0);
        elasticMatrix=[];
    end
    target = elasticMatrix;
    
elseif  flag ==  7 % atomic forces
    try
        [nothing, numStr] = unix('sed -n "6p" POSCAR');
        numIons= sum( str2num(numStr) );
        force_orig = callAWK('VASP_force.awk', 'OUTCAR', ['num=', num2str(numIons)]);
        if isempty(force_orig)
            USPEXmessage(1109,'',0);
            target = [];
        else
            target = force_orig;
        end
    catch
        target = [];
    end
elseif flag ==  8 % atomic positions in XDATA
    try
        [nothing, numStr] = unix('sed -n "6p" POSCAR');
        numIons= sum( str2num(numStr) );
        coord_orig = callAWK('VASP_xdatcar.awk', 'XDATCAR', ['num=', num2str(numIons)]);
        if isempty(coord_orig)
            USPEXmessage(1110,'',0);
            target = [];
        else
            target = coord_orig;
        end
    catch
        target = [];
    end
elseif flag ==  9 % lattice parameter in XDATA
    try
        lattice_orig = callAWK('VASP_xdatcar.awk', 'XDATCAR');
        if isempty(coord_orig)
            target = [];
        else
            target = lattice_orig(4:end, :);
        end
    catch
        target = [];
    end    
end

