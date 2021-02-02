function [potentialOffspring, template_for_random] = Random_Init_M400(angles_num)

potentialOffspring = [];

num_sec_str=fix(angles_num/7);
k = randi([1+fix(angles_num/49),num_sec_str],1,1);

l=zeros(1,k);

while l(end) <= 0
   for i=1:(size(l,2)-1)
      %r=randint(1,1,[1,fix(angles_num/k/2.5)])
      r=randi([1,fix(angles_num/k/3)],1,1);
      l(i)=randi([fix(angles_num/k)-r;fix(angles_num/k)+r],1,1);
   end

   l(end)=angles_num-sum(l);
end
%potentialOffspring = [potentialOffspring; [double(phi), double(psi)]]
l;
names_all=[];
All_res=dlmread('res');
for i=1:size(l,2)
   which_res=0;
   % Нам нужно, чтобы длина l, который мы хотим вытащить из pdb файла была меньше либо равен чем длина всего белка (which_res)
   while which_res < l(i)
      l(i)
      RandomIndex=randi([1, length(All_res)], 1,1)
      which_res = All_res(RandomIndex)
   end
   %which_res = randi([max(l(i),8), 49],1,1);
   folder=['/home/prachitskiy/DB/' num2str(which_res) '_res'];
   %folder=['./DB/' num2str(which_res) '_res'];
   all_pdb=dir(folder);
   try
       which_pdb=randi([3, size(all_pdb,1)], 1,1);
   catch ME
       causeException = MException('MATLAB:file', folder);
       ME = addCause(ME,causeException);
       rethrow(ME)
   end
   all_pdb(which_pdb);
   [names, angles] = Sec_from_stride([folder '/' all_pdb(which_pdb).name], l(i)); 
   names_all=[names_all names]
   %disp(names');
   potentialOffspring = [potentialOffspring; angles];
   template_for_random = all_pdb(which_pdb).name;
end

disp(names_all');
end
