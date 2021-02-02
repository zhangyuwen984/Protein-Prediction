function [coor, lat] = Read_QE_Structure()
%This rountine is to read structure from QE
%File: output (4.6/5.2)
%Last updated by Qiang Zhu (2013/10/03)


bohr = 0.529177; % Angstrom

[nothing, ans] = unix('grep CELL_P output |tail -1');
if ~isempty(findstr('alat',ans))  %The line must have alat keywords
%CELL_PARAMETERS (alat=  8.98730000)
%   0.999526127   0.000354354  -0.002205037
%  -0.235971977   0.956969070   0.000738780
%   0.148207110  -0.140193788   0.555468321   
   tmp = ans(23:32);
   scale_fac = str2num(tmp);
   %[a,LAT]=unix('awk -f QE_cell.awk output |tail -3');
   %lat = str2num(LAT);
   lat = callAWK('QE_cell.awk', '| tail -3', 'output');
else  
%crystal axes: (cart. coord. in units of a_0)
%          a(1) = (   1.000000   0.000000   0.000000 )
%          a(2) = (   0.000000   1.200794   0.000000 )
%          a(3) = (  -0.424832  -0.433458   0.820128 )
   [nothing, scaleStr] = unix('./getStuff output "lattice parameter" 6');
   scale_fac = str2num(scaleStr);
   %[a,LAT]=unix('awk -f QE_cell.awk output |tail -3');
   %lat = str2num(LAT);
   lat = callAWK('QE_cell.awk', '| tail -3', 'output');
end
lat = lat*bohr;
lat = lat*scale_fac;

[a,SCF]=unix('grep ATOMIC_POSITIONS output |wc -l');  %
scf = str2num(SCF);

%ATOMIC_POSITIONS (crystal)
%C        0.543835497   0.918532250   0.208306261
%C        0.596276457   0.404105264   0.975908592
%C        0.881126213   0.565692840   0.022620716
%C        0.450922873   0.080898506   0.853660863
%[a,COOR]=unix('awk -f QE_atom.awk output');
%COOR = str2num(COOR);
COOR = callAWK('QE_atom.awk','output');

numIons = size(COOR,1)/scf;
if scf > 1
   COOR(1:numIons*(scf-1),:)=[];
end
coor = COOR;
