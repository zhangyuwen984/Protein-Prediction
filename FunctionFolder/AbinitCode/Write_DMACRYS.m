function Write_DMACRYS(Ind_No)
%%%%%%%%%%%%%%%%%%%added on Mar 16th,2010%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global POP_STRUC
global ORG_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial files: .res, .punch, .axes, .pots, cutoff,fort22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = POP_STRUC.POPULATION(Ind_No).Step;
numMols = POP_STRUC.POPULATION(Ind_No).numMols;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
MtypeLIST = POP_STRUC.POPULATION(Ind_No).MtypeLIST;
typesAList = POP_STRUC.POPULATION(Ind_No).typesAList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%1st step: prepare .res, punch, need to read information of coords and lat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [nothing, nothing] = unix(['cat /dev/null > mol.res.punch']);
 [nothing, nothing] = unix(['cat /dev/null > mol.res']);
      fid = fopen('mol.res','wt');
      fid2 = fopen('mol.res.punch','wt');
      fidm = fopen('mol.axes','wt');
      fprintf(fidm,'MOLX ');
      for i=1:length(numMols)
           fprintf(fidm,'%2g', numMols(i));
      end
      fprintf(fidm,'\n');
      
      fprintf(fid,'TITL STURCTURE\n');
      fprintf(fid2,'! XIVSFM\n');     
      fprintf(fid2,'\n');
    if ~(length(POP_STRUC.POPULATION(Ind_No).LATTICE)==6)
         lat=latConverter(POP_STRUC.POPULATION(Ind_No).LATTICE);
    else
        lat=POP_STRUC.POPULATION(Ind_No).LATTICE;
    end
     lat(4:6)=lat(4:6)*180/pi;
     fprintf(fid,'CELL 1.0000 %8g %8g %8g %8g %8g %8g\n', lat);
     fprintf(fid,'LATT -1\n');  
     fprintf(fid,'SFAC  '); 
       for i=1:length(ORG_STRUC.atomType)
           fprintf(fid,'%4s', megaDoof(ORG_STRUC.atomType(i)));
       end
     fprintf(fid,'\n');
  
     num=zeros(1,length(typesAList));
     newCOORDS=[];
        for checkInd=  1: sum(numMols)
            newCOORDS= cat(1,newCOORDS,POP_STRUC.POPULATION(Ind_No).MOLECULES(checkInd).MOLCOORS);
        end
      newCOORDS1=newCOORDS/POP_STRUC.POPULATION(Ind_No).LATTICE;
someInd=0;
for i=1:length(MtypeLIST)
    for j=1:length(ORG_STRUC.STDMOL(MtypeLIST(i)).types)
           someInd = someInd+1;
            type = megaDoof(typesAList(someInd));
             for ind=1:length(ORG_STRUC.atomType)
                 if typesAList(someInd)==ORG_STRUC.atomType(ind)
                    num(ind)=num(ind)+1;
                    index = ind;
                 end
             end
            fprintf(fid,'%1s%1d %8g %9.5g %9.5g %9.5g\n',type,num(index),index,newCOORDS1(someInd,:) );
            if ~isempty(ORG_STRUC.STDMOL(MtypeLIST(i)).symbol)
                symbol = ORG_STRUC.STDMOL(MtypeLIST(i)).symbol(j);
            else
                symbol = '1';
            end
            label{someInd}=strcat(type,'_F',symbol,'_',num2str(num(index)),'____');
            if someInd==length(typesAList)
                  fprintf(fid2,'%3d %10s %10.5g %10.5g %10.5g Next   0 Limit  0\n ',someInd,label{someInd},newCOORDS(someInd,:) );
            else
                  fprintf(fid2,'%3d %10s %10.5g %10.5g %10.5g Next %4d Limit  0\n ',someInd,label{someInd},newCOORDS(someInd,:),someInd+1 );
            end
            if j==3
            fprintf(fidm,'X LINE %10s %10s 1 \n',label{someInd-2}, label{someInd-1} );
            fprintf(fidm,'Y PLANE %10s %10s 1 %10s 1\n',label{someInd-2}, label{someInd-1}, label{someInd} );
            end 
            fprintf(fid2,'0.0\n' );
            fprintf(fid2,'\n' );
      end   
end

            fprintf(fidm,'ENDS\n', numIons);
     fclose(fid);
     fclose(fid2);
     fclose(fidm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%2nd step: use neighcrys to generate .dmain and fort 20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    unix(['./neighcrys < fort.22']);
%if ORG_STRUC.constLattice
%          unix(['./neighcrys-vv < fort.22']);
%else
%    if(mod(step,2)==1)
%          unix(['./neighcrys-vv < fort.22']);
%    else
%          unix(['./neighcrys-pp < fort.22']);
%    end
%end
