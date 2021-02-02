function initialize_POP_STRUC_M300()

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{}, 'DoneOrder',{}, ...
                   'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'convex_hull',{}, 'fitness', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons', {},...
  'struc_entr',{}, 'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {},'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
                                                   'S_order',{},   'howCome',{},'JobID',{},'Folder',{}, 'typesAList',{},'chanAList',{},...
                                                'GB_COORDINATES',{}, 'GB_LATTICE',{}, 'GB_order',{},'GB_numIons',{},'GB_typesAList',{},...
                                   'Bulk_COORDINATES',{},'Bulk_LATTICE',{}, 'Bulk_numIons',{},'Bulk_typesAList',{}, 'Number',{}, 'symg',{});
POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{});


POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;
bulk_lat    =ORG_STRUC.bulk_lat;
bulk_pos    =ORG_STRUC.bulk_pos;
bulk_atyp   =ORG_STRUC.bulk_atyp';
bulk_numIons=ORG_STRUC.bulk_ntyp;
Lat         = latConverter(bulk_lat);
Lat(3)      = ORG_STRUC.thicknessS;
GB_numIons=ORG_STRUC.numIons;
item = 0;
for i = 1:length(GB_numIons)
    GB_atyp(1,item+1:item+GB_numIons(i))=ORG_STRUC.atomType(i);
    item = item + GB_numIons(i);
end

numIons = GB_numIons + bulk_numIons;

   for inPop_loop = 1:ORG_STRUC.initialPopSize
       POP_STRUC.POPULATION(inPop_loop).Bulk_LATTICE    =bulk_lat;
       POP_STRUC.POPULATION(inPop_loop).Bulk_COORDINATES=bulk_pos;
       POP_STRUC.POPULATION(inPop_loop).Bulk_typesAList =bulk_atyp;
       POP_STRUC.POPULATION(inPop_loop).Bulk_numIons    =bulk_numIons;
   end

    goodPop = 1;
    nsym1 = 0;

    while goodPop < ORG_STRUC.initialPopSize+0.5
%%%%%Applying Symmterization
      goodBad = 0;
      newSym = 0;
      while ~newSym
         newSym = 1;
         tmp = find(ORG_STRUC.nsym > 0);
         nsym = tmp(ceil(rand*length(tmp)));
         if nsym == nsym1
            newSym = 0;
         end
      end
      lat1 = latConverter(Lat);

      if sum(ORG_STRUC.splitInto)<2
         cd(['CalcFoldTemp']);
         [GB_candidate, GB_lat, errorS] = symope_GB(nsym, GB_numIons, lat1, ORG_STRUC.minDistMatrice);
         cd ..
      else
         splitInto = ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand));
         [GB_candidate, GB_lat, split, errorS] = splitBigCell_M300(lat1, splitInto, GB_numIons, nsym, ORG_STRUC.minDistMatrice);
      end
      if errorS == 0
         [lat,candidate, typesAList, chanAList] = makeGB(numIons, GB_lat, GB_candidate, GB_atyp, bulk_lat, bulk_pos, bulk_atyp, ORG_STRUC.vacuumSize(1));
         goodBad = distanceCheck(candidate, lat, numIons, ORG_STRUC.minDistMatrice);
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if goodBad
           nsym1 = nsym;
           POP_STRUC.POPULATION(goodPop).COORDINATES = candidate;
           POP_STRUC.POPULATION(goodPop).numIons = numIons;
           POP_STRUC.POPULATION(goodPop).LATTICE = lat;
           POP_STRUC.POPULATION(goodPop).typesAList = typesAList;
           POP_STRUC.POPULATION(goodPop).GB_LATTICE = GB_lat;
           POP_STRUC.POPULATION(goodPop).GB_COORDINATES = GB_candidate;
           POP_STRUC.POPULATION(goodPop).GB_numIons = GB_numIons;
           POP_STRUC.POPULATION(goodPop).GB_typesAList = GB_atyp;
           POP_STRUC.POPULATION(goodPop).chanAList = chanAList;
           POP_STRUC.POPULATION(goodPop).howCome = '  Random  ';

           if sum(ORG_STRUC.splitInto)>=2
              disp(['split into: ' num2str(split) '  subcells']);
           end

           disp(['Structure ' num2str(goodPop) ' with plane group ' num2str(nsym) ' generated randomly']);

%unix(['echo group_' num2str(nsym) ' >> SymmetryPOSCAR']);
%unix('echo 1.0 >> SymmetryPOSCAR');
%unix(['echo ' num2str(lat(1,:),'%12.3f') ' >> SymmetryPOSCAR']);
%unix(['echo ' num2str(lat(2,:),'%12.3f') ' >> SymmetryPOSCAR']);
%unix(['echo ' num2str(lat(3,:),'%12.3f') ' >> SymmetryPOSCAR']);
%unix(['echo ' num2str(numIons) ' >> SymmetryPOSCAR']);
%unix('echo Direct >> SymmetryPOSCAR');
%for i = 1 : sum(numIons)
%unix(['echo ' num2str(candidate(i,:),'%12.3f') ' >> SymmetryPOSCAR']);
%end

           goodPop = goodPop+1;

        end
    end

Start_POP_M300();


