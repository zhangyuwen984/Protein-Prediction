function kuerzel = valence(Ordnungszahl) % number of valence electrons

%nVal = cellstr(['1';'2';'1';'2';'3';'4';'3';'2';'1';'1';'1';'2';'3';'4';'3';'2';'1';'1';'1';'2';...
%   '3';'4';'5';'6';'5';'3';'3';'3';'2';'2';'3';'4';'3';'2';'1';'1';'1';'2';'3';'4';...
%   '5';'6';'5';'3';'3';'3';'2';'2';'3';'4';'3';'2';'1';'1';'1';'2';'3';'3';'3';'3';...
%   '3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';...
%   '3';'4';'3';'2';'1';'3';'1';'2';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';'3';...
%   '3';'3';'3';'3';'3']);


nVal = cellstr([
   '1.0';'0.5';'1.0';'2.0';'3.0';'4.0';'3.0';'2.0';'1.0';'0.5';'1.0';'2.0';'3.0';'4.0';'3.0';'2.0';'1.0';'0.5';'1.0';'2.0';...
   '3.0';'4.0';'4.0';'3.0';'4.0';'3.0';'3.0';'2.0';'2.0';'2.0';'3.0';'4.0';'3.0';'2.0';'1.0';'0.5';'1.0';'2.0';'3.0';'4.0';...
   '5.0';'4.0';'4.0';'4.0';'4.0';'4.0';'1.0';'2.0';'3.0';'4.0';'3.0';'2.0';'1.0';'0.5';'1.0';'2.0';'3.0';'4.0';'3.0';'3.0';...
   '3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'3.0';'4.0';'5.0';'4.0';'4.0';'4.0';'4.0';'4.0';'1.0';'2.0';...
   '3.0';'4.0';'3.0';'2.0';'1.0';'0.5';'1.0';'2.0';'3.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';'4.0';...
   '4.0';'4.0';'4.0';'4.0';'2.0';'2.0';'2.0';'2.0';'2.0';'2.0']);


kuerzel = nVal{Ordnungszahl};
