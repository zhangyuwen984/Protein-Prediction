function uff1 = ReviseMOL(uff, atomType)

uff1 = uff;
tmp = find(uff.bound==1);
if length(tmp) ~=2
      disp('For polymers, there must exist 2 active atoms in MOL_1 ...')
      disp('Please go back and recheck the MOL_1 file...')
      quit;
else
      format = uff.format;
      coords = uff.molecule;
      bound  = uff.bound;
      symbol = uff.symbol;
      types  = uff.types;
      charge = uff.charge;

    if tmp(1) ~= 1
      disp('We need to resort the atomic order for polymers...');

      for i = 1:length(types)
          radii(i) = str2num(covalentRadius(atomType(types(i))));
      end
      %Step2: find atomic pair
       Pair = find_pair(coords, radii);
      %Step3, resort the atom, by placing 1st bound at the 1st place
       [coords, symbol, charge, types, bound, Pair] = ...
          resortXYZ(coords, symbol, charge, types, bound, Pair, tmp(1), radii);
      %Step4: find the format
       format = findformat(Pair);
    end
    z=NEW_coord2Zmatrix(coords,format);
    z(1,:)=[0 0 0];
    z(2,2)=0;
    z(2,3)=0;
    z(3,3)=0;
    coord1=[];
    [coord1, Len] = MergeZmatrix(z, format, bound, atomType(types));
    uff1.molecule = coord1;
    uff1.format   = format;
    uff1.symbol   = symbol;
    uff1.charge   = charge;
    uff1.types    = types;
    uff1.bound    = bound;
    uff1.length   = Len;
    
end

