function [lattice,coordinates,at,type]=supercell(lattice, coordinates, atomtype, ntyp, cell)
         
       global ORG_STRUC 
         cor=[];
         atomType=[];
         LAT=latConverter(lattice);
         if cell(1)<1
            coor1=coordinates(:,1)*LAT(1);
            standard=LAT(1)*cell(1);
            k=0;
            for i=1:size(coor1,1)
                if coor1(i)<=standard
                   k=k+1;
                   cor(k,:)=coordinates(i,:);
                   cor(k,1)=coor1(i)/standard;
                   atomType(k)=atomtype(i);
                end
            end
         else
            for i=1:cell(1)
              POS=coordinates;
              POS(:,1)=coordinates(:,1)+i-1;
              cor=cat(1,cor,POS);
              a=size(atomtype);
              if(a(1)<a(2))
                 atomtype=atomtype';
              end
              atomType=cat(1,atomType,atomtype);
            end
            cor(:,1)=cor(:,1)/cell(1);
          end
         cor1=cor;
         atomType1=atomType;
         atomType2=[];
         cor2=[];
          if cell(2)<1
            coor2=cor1(:,2)*LAT(2);
            standard=LAT(2)*cell(2);
            k=0;
            for i=1:size(coor2,1)
                if coor2(i)<=standard
                   k=k+1;
                   cor2(k,:)=cor1(i,:);
                   cor2(k,2)=coor2(i)/standard;
                   atomType2(k)=atomType1(i);
                end
            end
         else
            for j=1:cell(2)
                POS=cor1;
                POS(:,2)=cor1(:,2)+j-1;
                cor2=cat(1,cor2,POS);
                a=size(atomType1);
                if(a(1)<a(2))
                   atomType1=atomType1';
                end
                atomType2=cat(1,atomType2,atomType1);
            end
            cor2(:,2)=cor2(:,2)/cell(2);
         end
         cor3=cor2;
         atomType3=atomType2;
         cor4=[];
         atomType4=[];
         if cell(3)<1
            coor3=cor3(:,3)*LAT(3);
            standard=LAT(3)*cell(3);
            k=0;
            for i=1:size(coor3,1)
                if coor3(i)<=standard
                   k=k+1;
                   cor4(k,:)=cor3(i,:);
                   cor4(k,3)=coor3(i)/standard;
                   atomType4(k)=atomType3(i);
                end
            end
         else
            for k=1:cell(3)
              POS=cor3;
              POS(:,3)=cor3(:,3)+k-1;
              cor4=cat(1,cor4,POS);
              a=size(atomType3);
              if(a(1)<a(2))
                   atomType3=atomType3';
              end
              atomType4=cat(1,atomType4,atomType3);
            end
            cor4(:,3)=cor4(:,3)/cell(3);
         end
         cord=cor4;
         at=atomType4;
         lat=latConverter(lattice);
         lat(1)=lat(1)*cell(1);
         lat(2)=lat(2)*cell(2);
         lat(3)=lat(3)*cell(3);
         lattice=latConverter(lat);
         coordinates=cord;
         c=size(atomType4);
         if c(1)>c(2)
            atomType4=atomType4';
         end
         for i=1:size(ORG_STRUC.atomType,2)
             type(i)= sum(size(find(atomType4==ORG_STRUC.atomType(i)),2));
         end
