function [lat1,Surcandidate,vol,SuratomType,surnumIons,flag]=reduce_Surface(lattice,coordinates,Ind)
%The new idea is implemented to improve the efficiency. 
%Sometimes, the surface atoms will get far away from the bulk, 
%and then some extreme horriable structures appear. 
%To avoid this, we firstly fix all the bulk atoms in the first step. 
%If the surface atoms are not bonded but go to the vacuum, just delete them.   
     global ORG_STRUC 
     global POP_STRUC
     flag=0;
     Step=POP_STRUC.POPULATION(Ind).Step;
     MaxStep = length([ORG_STRUC.abinitioCode]);
     chanAList=POP_STRUC.POPULATION(Ind).chanAList;
     thicknessS = ORG_STRUC.thicknessS;
     atomtype=POP_STRUC.POPULATION(Ind).typesAList; 

     lat=latConverter(lattice);
     % find the Volume of small cell for fingerprint 
     maxz=max(coordinates(:,3)*lat(3));
     lat_print=[lat(1) lat(2) maxz  lat(4) lat(5) lat(6)];
     lattice_print=latConverter(lat_print);
     vol=det(lattice_print);
     Netcandidate=[];  %%%to store the add atoms in the network
     NetatomType=[];   %%%
%%%%%%%%%%%%%%%%%%
     if Step < (MaxStep + 0.5)
        k=0;
        [isAllConnected, isConnectedToSubstrate] = ...
        surface_connectivity_check( lattice, coordinates, atomtype, chanAList);
        if ~isAllConnected
           disp(['Structure ' num2str(Ind) ' not AllConected at Step ' num2str(Step)]);
           flag = 1;
        end
%%%%%%%just take the surface atoms in the network    
        for i=1:size(chanAList,2)
           if chanAList(i) & (isConnectedToSubstrate(i))
              k=k+1;
              Netcandidate(k,:)=coordinates(i,:);
              NetatomType(k)=atomtype(i);
           end
        end
      else
        k=0;
        for i=1:size(chanAList,2)
           if chanAList(i)
              k=k+1;
              Netcandidate(k,:)=coordinates(i,:);
              NetatomType(k)=atomtype(i);
           end
        end
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Surcandidate=[];
        SuratomType=[];
     if isempty(Netcandidate)        
        lat1=lat;
        lat1(3) = thicknessS;
        lat1=latConverter(lat1);
     else
        Netcandidate(:,3)=Netcandidate(:,3)*lat(3);
        dz=ORG_STRUC.bulk_lat(3,3);  %here we use bulk length as the lower bound
        if (Step < MaxStep+0.5)
            item=0;
            for i=1:k
               if (Netcandidate(i,3) < dz+thicknessS) & (Netcandidate(i,3)>dz-1)
                  flag=1;
                  item=item+1;
                  Surcandidate(item,:)=Netcandidate(i,:);
                  SuratomType(item)=NetatomType(i);
               else
                  disp('delete the atoms going outside the surface');
               end                    
            end
        else
            Surcandidate=Netcandidate;
            SuratomType=NetatomType;
        end

        lat(3)=thicknessS;
        for i=1:size(Surcandidate,1)
          Surcandidate(i,3)=Surcandidate(i,3)-dz;
        end
        lat=[lat(1) lat(2) lat(3)  lat(4) lat(5) lat(6)];
        lattice=latConverter(lat);
        for i=1:size(Surcandidate,1)
          Surcandidate(i,3)=Surcandidate(i,3)/lat(3);
        end
        lat1=lattice;
    end

    surnumIons=[];
    for i=1:length(ORG_STRUC.atomType)
         surnumIons(1,i)=sum( SuratomType==ORG_STRUC.atomType(i));    
    end
