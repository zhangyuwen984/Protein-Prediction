function writeCelFiles(Coordinates, Lattice, numIons, atomType, bodyCount, whereToWrite)

% USPEX Version 6.6.7
% Change: moved files to separate folder

global POP_STRUC

toWriteFile = [whereToWrite 'AuxiliaryFiles/OUTFile' num2str(bodyCount) '.CEL'];
[nothing, nothing] = unix(['cat /dev/null > ' toWriteFile]);

if size(Lattice,1) == 3
  Lattice = latConverter(Lattice);
end
if size(Lattice,1) == 1
  Lattice = Lattice';
end
Lattice(4:6) = Lattice(4:6)*180/pi;
[nothing, nothing] = unix(['echo cell ' num2str(Lattice') ' >> ' toWriteFile]);

[nothing, nothing] = unix(['echo natom ' num2str(sum(numIons)) ' >> ' toWriteFile]);

coordLoop = 1;
for i = 1 : length(numIons)
    sameAtoms = 1;
    atNum = num2str(atomType(i));
    atShort =  megaDoof(atomType(i));
    
    for j = 1 : numIons(i)
      [nothing, nothing] = unix(['echo ' atShort num2str(sameAtoms) ' ' atNum ' ' num2str(Coordinates(coordLoop,:)) ' >> ' toWriteFile]);
      sameAtoms = sameAtoms + 1;
      coordLoop = coordLoop + 1;
    end
end
[nothing, nothing] = unix(['echo rgnr 1 >> ' toWriteFile]);
