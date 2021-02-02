function writeCalcFilesTINKER_MOL(Ind_No)

global POP_STRUC
global ORG_STRUC

lattice    = POP_STRUC.POPULATION(Ind_No).LATTICE;
numMols    = POP_STRUC.POPULATION(Ind_No).numMols;
typesAList = POP_STRUC.POPULATION(Ind_No).typesAList;
MtypeLIST  = POP_STRUC.POPULATION(Ind_No).MtypeLIST;
MOLECULES  = POP_STRUC.POPULATION(Ind_No).MOLECULES;
STDMOL     = ORG_STRUC.STDMOL;
newCOORDS = [];
atomsMol  = [];
for checkInd = 1 : sum(numMols)
    % Counting number of atoms in each molecule:
    atomsMol  = [atomsMol, length(MOLECULES(checkInd).MOLCOORS)];
    newCOORDS = [newCOORDS; MOLECULES(checkInd).MOLCOORS];
end
newCOORDS = newCOORDS / lattice;

%fp = fopen('input.key','a+');
%fprintf(fp, 'spacegroup  P1\n');
%Lattice = latConverter(lattice);
%Lattice(4:6) = Lattice(4:6)*180/pi;
%fprintf(fp,'a-axis      %8.4f\nb-axis      %8.4f\nc-axis      %8.4f\nalpha       %8.4f\nbeta        %8.4f\ngamma       %8.4f\n', Lattice);
%fclose(fp);

%unix(['rm -f input.fra input.xyz* ']);
fp = fopen('input.fra','w');
fprintf(fp,'%4i atoms\n',length(typesAList) );
Lattice = latConverter(lattice);
Lattice(4:6) = Lattice(4:6)*180/pi;
fprintf(fp,'%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', Lattice);

lastAtom = 0;
cov1=zeros(length(typesAList), 4);
for someInd = 1 : length(typesAList)
    [i,j]=GetID(someInd, MtypeLIST, typesAList);
    kind=STDMOL(MtypeLIST(i)).charge(j);
    cov1(someInd, :)=STDMOL(MtypeLIST(i)).CN(j,1:4);
    cov1(someInd, 1) = cov1(someInd, 1) + lastAtom;
    
    if cov1(someInd, 2)==0
    else
        cov1(someInd, 2)=lastAtom + cov1(someInd, 2);
    end
    
    if cov1(someInd, 3)==0
    else
        cov1(someInd, 3)=lastAtom + cov1(someInd, 3);
    end
    
    if cov1(someInd, 4)==0
    else
        cov1(someInd, 4)=lastAtom + cov1(someInd, 4);
    end
    
    totalAtoms = 0;
    for checkInd = 1 : sum(numMols)
        totalAtoms = totalAtoms + atomsMol(checkInd);
        if someInd == totalAtoms
            lastAtom = someInd;
        end
    end
    %   nat=STDMOL(MtypeLIST(i)).nat(j);
    label = [megaDoof(typesAList(someInd))];
    
    if ~isempty(STDMOL(MtypeLIST(i)).symbol)
        symbol=STDMOL(MtypeLIST(i)).symbol(j);
        label = [label symbol]';
    end
    fprintf(fp,'%6i %4s %8.4f %8.4f %8.4f  %8i %8i %8i %8i %8i\n', someInd, label, newCOORDS(someInd,:), kind, cov1(someInd, :));
end
fclose(fp);

