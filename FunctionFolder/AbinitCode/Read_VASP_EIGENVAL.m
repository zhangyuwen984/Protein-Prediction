function [eGap] = read_EIGENVAL

% this function tries to determine the bandgap from the EIGANVAL file
% should be more exact than from DOSCAR for gapped materials
% Line 6: number of electrons, N of K-poins, N of bands
% every k-point: k-point coordinate, weight; then band index, eigenvalue 
% (for spin-polarized systems - 2 columns of eigenvalues for spin up and down respectively).

try
 handle = fopen('EIGENVAL');
 for i = 1 : 6
   tmp = fgetl(handle);
 end
 l6 = str2num(tmp);
 eN = l6(1);
 filledBand = floor(eN/2); % 2 electrons per band
 nK = l6(2); % number of k-points
 nB = l6(3); % number of bands
 for i = 1 : nK
   tmp = fgetl(handle); % empty line
   tmp = fgetl(handle); % k point coords and weight
   bands = fscanf(handle, '%g %g', [2 nB]); % do something about spin-up and down!!!!!
   bands = bands';
   if i == 1
     gapStart = bands(filledBand, 2);
     gapEnd = bands(filledBand+1, 2);
   else
     if gapStart < bands(filledBand, 2);
       gapStart = bands(filledBand, 2);
     end
     if gapEnd > bands(filledBand+1, 2);
       gapEnd = bands(filledBand+1, 2);
     end
   end
   tmp = fgetl(handle); % EOL char
 end
 fclose(handle);

 eGap = gapEnd - gapStart;
 if eGap < 0
   eGap = 0; % for metals
 end
catch
 eGap = 0;
end