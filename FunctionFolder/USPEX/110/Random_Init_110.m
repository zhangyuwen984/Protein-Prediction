function [Molecules, lat, MtypeLIST, typesAList, numIons] = Random_Init_110(Ind_No, numMols)

% implemented - USPEX Version 9.3.7
global ORG_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newSym = 1;
badSymmetry = 0;
goodBad = 0;
Molecules = struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{});

while ~goodBad

  [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
   if ORG_STRUC.constLattice    % lat1 = lattice
      lat1 = ORG_STRUC.lattice;
   else                         % lat1 = lattice volume
      lat1 = ORG_STRUC.latVolume/ORG_STRUC.STDMOL(1).length;
   end
   if badSymmetry > 15
      badSymmetry = 0;
      newSym = 1;      % change the symmetry group if can't generate the crystal
   end

   badSymmetry = badSymmetry + 1;
   if newSym
      tmp = find(ORG_STRUC.nsym > 0);
      nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
      newSym = 0;
   end
   CenterMat = ORG_STRUC.CenterminDistMatrice;

    if nsym > 1  % if nsym == 1 - just keep the random candidate that we generated
        cd([ORG_STRUC.homePath '/CalcFoldTemp'])
        %[candidate, lat, errorS] = symope_2D(nsym, numMols, lat1, CenterMat, 0);
        [candidate, lat, errorS] = symope_2D(nsym, numMols, lat1, CenterMat, 0);
        lat(3,3) = ORG_STRUC.STDMOL(1).length;
        cd(ORG_STRUC.homePath)
    end

    if errorS==0
      for item=1:30
        Molecules = [];
        %count = 1;
        for ind= 1: sum(numMols)
            format   = ORG_STRUC.STDMOL(MtypeLIST(ind)).format;
            MOLCOORS = ORG_STRUC.STDMOL(MtypeLIST(ind)).molecule;
            MOLCOORS = bsxfun(@minus, MOLCOORS, mean(MOLCOORS));
            if rand > 0.5
               MOLCOORS = -1*MOLCOORS;
            end
            R_axis = [0 0 1];
            angle=rand*pi;
            MOLCOORS = Rotate_rigid_body([0 0 0], R_axis, MOLCOORS, angle);
            MOLCOORS = bsxfun(@plus, MOLCOORS, (candidate(ind,:)+[0 0 rand*0.5])*lat);
            Molecules(ind).MOLCOORS= MOLCOORS;
            Molecules(ind).ZMATRIX = real(NEW_coord2Zmatrix(MOLCOORS,format));
           % count = count + 1;
           % if innerInder > 2
           %    angle=rand*pi;
           %    MOLCOORS = Rotate_rbody(mean(MOLCOORS), R_axis, MOLCOORS, angle);
           %    MOLCOORS = bsxfun(@minus, MOLCOORS, mean(MOLCOORS));
           %    MOLCOORS = bsxfun(@plus, MOLCOORS, (candidate(innerInder,:)+ [0 0 0.5])*lat);
           %    Molecules(count).MOLCOORS= MOLCOORS;
           %    Molecules(count).ZMATRIX = real(NEW_coord2Zmatrix(MOLCOORS,format));
           %    count = count + 1;
           % end
        end
           goodBad = newMolCheck(Molecules,lat, MtypeLIST, ORG_STRUC.minDistMatrice);
        if goodBad
           disp(['Structure ' num2str(Ind_No) ' bulit with the plane group ' num2str(nsym)]);
           break
        end
      end
     end
  end
