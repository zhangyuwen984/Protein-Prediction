function Dummy = QuickStart(INPUT)
%To set all the fields empty
%Octave sometimes fails to keep the field for the struc
%After we declare it, we must initialize it
%Qiang Zhu(2013/01/14)
Names = fieldnames(INPUT);
Dummy = struct();
Num = size(Names, 1);
for i = 1:Num
   Dummy = setfield(Dummy, Names{i}, '');
end
