function Write_GULP_MOL(Ind_No, Write_Conn)

global POP_STRUC
global ORG_STRUC

lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
numMols = POP_STRUC.POPULATION(Ind_No).numMols;
typesAList= POP_STRUC.POPULATION(Ind_No).typesAList;
MtypeLIST= POP_STRUC.POPULATION(Ind_No).MtypeLIST;
MOLECULES= POP_STRUC.POPULATION(Ind_No).MOLECULES;
STDMOL   = ORG_STRUC.STDMOL;
newCOORDS = [];
for checkInd = 1 : sum(numMols)
    newCOORDS = cat(1,newCOORDS, MOLECULES(checkInd).MOLCOORS);
end
    newCOORDS = newCOORDS / lattice;

fp = fopen('input','a+');
fprintf(fp, 'cell\n');
Lattice = latConverter(lattice);
Lattice(4:6) = Lattice(4:6)*180/pi;
fprintf(fp,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', Lattice);
fprintf(fp,'fractional\n');
for someInd = 1 : length(typesAList)
   [i,j]=GetID(someInd, MtypeLIST, typesAList);
   charge=STDMOL(MtypeLIST(i)).charge(j);
   label = [megaDoof(typesAList(someInd))];
   if ~isempty(STDMOL(MtypeLIST(i)).symbol)
       symbol=STDMOL(MtypeLIST(i)).symbol(j);
       label = [label '_' symbol]';
   end
   fprintf(fp,'%4s %8.4f %8.4f %8.4f core %8.6f\n', label, newCOORDS(someInd,:), charge);
end

%%%%%%%%%%Only for test
%if ORG_STRUC.varcomp == 1
%count = 0;
%   for i = 1 : length(MtypeLIST)
%       if MtypeLIST(i) == 1  %C-H bond in CH4
%          fprintf(fp, 'connect %4d %4d single 0 0 0\n', count+1, count+2);
%          fprintf(fp, 'connect %4d %4d single 0 0 0\n', count+1, count+3);
%          fprintf(fp, 'connect %4d %4d single 0 0 0\n', count+1, count+4);
%          fprintf(fp, 'connect %4d %4d single 0 0 0\n', count+1, count+5);
%       else  % H-H bond in H2
%          fprintf(fp, 'connect %4d %4d single 0 0 0\n', count+1, count+2);
%       end
%          count = count + length(STDMOL(MtypeLIST(i)).types);
%   end
%end

if Write_Conn
   count = 0;
   for i = 1 : length(MtypeLIST)
       for j = 1:length(STDMOL(MtypeLIST(i)).types)   %C-H bond in CH4
           CN = STDMOL(MtypeLIST(i)).CN(j,:); %neighbor list
           for k = 1:length(CN)-1  %the last CN element is the coordination number
               if (CN(k) > 0) & (CN(k) > j)
                  fprintf(fp, 'connect %4d %4d \n', count+j, count+CN(k));
               end
           end
       end
       count = count + length(STDMOL(MtypeLIST(i)).types);
   end
end

fclose(fp);
