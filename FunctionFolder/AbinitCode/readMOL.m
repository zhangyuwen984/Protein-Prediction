function readMOL(Ind_No, code)

%This rountine is to construct MOLECULE data structure for three cases
%0, The output is VASP POSCAR which counts element by element (For reading Seeds)
%1, The output is VASP POSCAR which counts element by element (For general case)
%2, The output is Zmatrix (for siesta)
%3, The output is xyz format which counts molecule by molecule (For DMACRYS and perhaps GULP in future)
%last updated by Qiang Zhu (2013/10/05)
global POP_STRUC
global ORG_STRUC


oldPOP = POP_STRUC.POPULATION(Ind_No);  % store the old structure info 
 
MtypeLIST = POP_STRUC.POPULATION(Ind_No).MtypeLIST;
numMols = POP_STRUC.POPULATION(Ind_No).numMols;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
COORDINATES= POP_STRUC.POPULATION(Ind_No).COORDINATES;
typesAList = POP_STRUC.POPULATION(Ind_No).typesAList;

if code==6
   maxdist = 0.05;
else
   maxdist = 0.6;
end

if code==6 %DMACRYS
   code=3;
elseif code==12 %TINKER
   code=3;
elseif code>3
   code=1;
end

%GET ZMATRIX, MOLCOORS, and COORDINATES
if code == 2  %FOR SIESTA, we read Zmatrix
     newCoords = COORDINATES;
     newCoords1 = zeros(0,3);
     for ind = 1: sum(numMols)
         nowMOL = MtypeLIST(ind);
         format = ORG_STRUC.STDMOL(nowMOL).format;
         amount = length(ORG_STRUC.STDMOL(nowMOL).types); %todo
         ZMATRIX= newCoords(1:amount,:);
         MOLCOORS = NEW_ZMATRIXCOORD(ZMATRIX, format);
         POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).ZMATRIX  =ZMATRIX;
         POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCOORS =MOLCOORS;
         newCoords(1:amount,:)=[];
         newCoords1 = cat(1,newCoords1, MOLCOORS);
         MOLCOORS = [];
         ZMATRIX  = [];
         if length(ORG_STRUC.STDMOL(MtypeLIST(ind)).types) ==1
            POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCENTER = MOLCOORS;
         else
            POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCENTER = mean(MOLCOORS);
         end
     end
     for m = 1: length(ORG_STRUC.atomType)
         s = find(typesAList == ORG_STRUC.atomType(m));
         newCoords1 = cat(1,newCoords1, newCoords(s,:));
     end
     POP_STRUC.POPULATION(Ind_No).COORDINATES = newCoords1;

else %For the other code, we read coodinates
   if code ==3  %For DMACRYS, TINKER and GULP
       newCoords = COORDINATES;
       newCoords1 = zeros(0,3);
       for m = 1: length(ORG_STRUC.atomType)
           s = find(typesAList == ORG_STRUC.atomType(m));
           newCoords1 = cat(1,newCoords1, newCoords(s,:));
       end
       POP_STRUC.POPULATION(Ind_No).COORDINATES = newCoords1;
   else %for VASP and general cases
      item = zeros(1,length(ORG_STRUC.atomType));
      newCoords = zeros(0,3);
      for mol = 1: length(MtypeLIST)
          for atom = 1: length(ORG_STRUC.STDMOL(MtypeLIST(mol)).types)
              ind = ORG_STRUC.STDMOL(MtypeLIST(mol)).types(atom);
              item(ind)=item(ind)+1;
              if ind > 1
                 s = sum(numIons(1:ind-1)) + item(ind);
              else
                 s = item(1);
              end
              newCoords = cat(1,newCoords, COORDINATES(s,:));
          end
      end
   end

   % To Update MOLECULES
   for ind = 1: sum(numMols)
       MOLCOORS = [];
       amount = length(ORG_STRUC.STDMOL(MtypeLIST(ind)).types); %todo
       MOLCOORS= newCoords(1:amount,:)*lattice;
       format = ORG_STRUC.STDMOL(MtypeLIST(ind)).format;
       MOLCOORS= Movecoor(MOLCOORS, lattice, format);
       POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCOORS= MOLCOORS;
       POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).ZMATRIX = real(NEW_coord2Zmatrix(MOLCOORS, format));
       if length(ORG_STRUC.STDMOL(MtypeLIST(ind)).types) ==1
          POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCENTER = MOLCOORS;
       else
          POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).MOLCENTER = mean(MOLCOORS);
       end
       newCoords(1:amount,:)=[];
       MOLCOORS = [];
   end
end

if ~isempty(POP_STRUC.POPULATION(Ind_No).Step) %except Seeds
   if ORG_STRUC.checkMolecules
      for ind = 1: sum(numMols)
          zmatrix = POP_STRUC.POPULATION(Ind_No).MOLECULES(ind).ZMATRIX;
          zmatrix_STD = ORG_STRUC.STDMOL(MtypeLIST(ind)).ZMATRIX;
          decompose = 0;
          if size(zmatrix,1) > 1
             for loop = 2:size(zmatrix,1)
                 if abs(zmatrix(loop,1)-zmatrix_STD(loop,1))>maxdist %bond length should not be too big!        
                    decompose = 1;
          	  break;
                 end
             end
          end

          Poly=newMolCheck(POP_STRUC.POPULATION(Ind_No).MOLECULES, lattice,MtypeLIST, ORG_STRUC.minDistMatrice);

          if (decompose == 1) | (Poly==0)
             POP_STRUC.POPULATION(Ind_No).Enthalpies(end)=100000;  %always reset the unsatisfied structure to 100000
             disp('MOLECULE FORM is not kept, see Warnings for more details ...');
             
             molePOSCARStr = createPOSCARStr(oldPOP, POP_STRUC.POPULATION(Ind_No), ORG_STRUC.atomType);
             USPEXmessage(0, molePOSCARStr, 1);
             if POP_STRUC.POPULATION(Ind_No).Step>=2
             	POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
             end
             break;
          end
      end
   end
end


function newCoor = Movecoor(coords,lat,format)

if(size(coords,1)>1)
   X1 = [-2:2];
   Y1 = [-2:2];
   Z1 = [-2:2];
   [X2,Y2,Z2] = meshgrid(X1,Y1,Z1);
   X3=reshape(X2,1,[]);
   Y3=reshape(Y2,1,[]);
   Z3=reshape(Z2,1,[]);
   Matrix = [X3; Y3; Z3]';
   for ind = 2:size(coords,1);
      coor1 = coords(format(ind,1),:);
      coor2 = coords(ind,:);
      tmp   = bsxfun(@plus, Matrix*lat, coor2);
      [min_dist,ID] = min(pdist2(tmp, coor1));
      coords(ind,:) = tmp(ID,:);
   end
end
newCoor = coords;


function molePOSCARStr = createPOSCARStr(oldPOP, POP, atomType )

atomTypStr = [];
for i = 1:length(atomType)
  atomTypStr = [atomTypStr, ' ', megaDoof( atomType(i) ) ' '];
end



% after Relaxation
molePOSCARStr= ['MOLECULE FORM is not kept in this POSCAR\n', '1.00\n'];

for i = 1:3
  molePOSCARStr= [ molePOSCARStr, num2str(POP.LATTICE(i,:)) '\n'];
end

molePOSCARStr= [ molePOSCARStr, atomTypStr, '\n'];
molePOSCARStr= [ molePOSCARStr, num2str(POP.numIons),'\n', 'Direct' '\n'];
for i = 1:sum( POP.numIons )
   molePOSCARStr= [ molePOSCARStr, num2str(POP.COORDINATES(i,:)) '\n' ];
end


% The original Structure 
molePOSCARStr= [molePOSCARStr, '\n\nStruture before relaxation, for checking\n', '1.00\n'];

for i = 1:3
  molePOSCARStr= [ molePOSCARStr, num2str(oldPOP.LATTICE(i,:)) '\n'];
end

molePOSCARStr= [ molePOSCARStr, atomTypStr, '\n'];
molePOSCARStr= [ molePOSCARStr, num2str(oldPOP.numIons),'\n', 'Direct' '\n'];
for i = 1:sum( oldPOP.numIons )
   molePOSCARStr= [ molePOSCARStr, num2str(oldPOP.COORDINATES(i,:)) '\n' ];
end

