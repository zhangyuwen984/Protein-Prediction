function [eGap, fermiDOS] = read_DOSCAR(nAtoms)

% this function tries to determine the bandgap from the DOSCAR file
% fermi energy - 4th value in the 6th string 

handle = fopen('DOSCAR');
for i = 1 : 6
  tmp = fgetl(handle);
end
line6 = str2num(tmp);
eFermi = line6(4);
DOS = fscanf(handle, '%g %g %g', [3 inf]);
DOS = DOS';
fclose(handle);

for eF = 1 : size(DOS,1)
  if DOS(eF,1) > eFermi
    break;
  end
end
% integrated DOS at fermi level, linear approximation:
fermiDOS = DOS(eF,3) - (DOS(eF,3) - DOS(eF-1,3))*(DOS(eF,1) - eFermi)/(DOS(eF,1) - DOS(eF-1,1)); 
k = (DOS(eF,3) - DOS(eF-1,3))/(DOS(eF,1) - DOS(eF-1,1)); % 1st derivative, linear approximation
k_next1 = (DOS(eF+1,3) - DOS(eF,3))/(DOS(eF+1,1) - DOS(eF,1));
k_next2 = (DOS(eF+2,3) - DOS(eF+1,3))/(DOS(eF+2,1) - DOS(eF+1,1));

if ((DOS(eF+1,3) - fermiDOS) > 0.01*nAtoms) | ((k_next1 > k) & (k_next2 > k))
  eGap = 0; % we got a metal
else
  k = (DOS(eF,3) - DOS(eF-1,3))/(DOS(eF,1) - DOS(eF-1,1)); % 1st derivative, linear approximation
  for i = eF + 1 : size(DOS,1)
    k1 = (DOS(i,3) - DOS(i-1,3))/(DOS(i,1) - DOS(i-1,1)); % 1st derivative, linear approximation
    if k1 > k
      break;
    end
  end
  eGap = DOS(i-1,1) - eFermi;
end

fermiDOS = (DOS(eF,2) - (DOS(eF,2) - DOS(eF-1,2))*(DOS(eF,1) - eFermi)/(DOS(eF,1) - DOS(eF-1,1)))/nAtoms;