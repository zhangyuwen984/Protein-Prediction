function Finish()
%add a complete sign, generally work for every type of calculation
%Last updated by Qiang Zhu (2013/10/17)

global POP_STRUC
global ORG_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

fprintf(fp,'This calculation runs %4d structural relaxations\n', POP_STRUC.bodyCount);

if isfield(POP_STRUC, 'convex_hull')
   if ~isempty(POP_STRUC.convex_hull)
      fprintf(fp,'The stable compositions are: \n');
      for i=1:size(POP_STRUC.convex_hull,1)
         fprintf(fp,'%4d %4d %8.4f\n', POP_STRUC.convex_hull(i,1:3));
      end
   else
      plotE(POP_STRUC.resFolder);
   end
else
   plotE(POP_STRUC.resFolder);
end

fprintf(fp,'\n', datestr(now));
fprintf(fp,'            Job Finished at       %30s\n', datestr(now));
POP_STRUC.generation = ORG_STRUC.numGenerations + 1;
safesave('Current_POP.mat', POP_STRUC)
disp(' ')
[nothing, nothing] = unix('rm still_running');
