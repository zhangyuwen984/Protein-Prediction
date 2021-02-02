function potentialOffspring = Random_Init_M400(angles_num)
% $Rev: 690 $
% $Author: maxim $
% $Date: 2014-10-28 07:25:19 +0400 (Tue, 28 Oct 2014) $

potentialOffspring = [];

num_sec_str=fix(angles_num/7);
k = randi([1,num_sec_str],1,1);

l=zeros(1,k);

for i=1:(size(l,2)-1)
   %r=randint(1,1,[1,fix(angles_num/k/2.5)])
   r=randi([1,fix(angles_num/k/3)],1,1);
   l(i)=randint(1,1,[fix(angles_num/k)-r;fix(angles_num/k)+r]);
end

l(end)=angles_num-sum(l);

for i=1:size(l,2)
    which_sec_str=randi([1,7],1,1);
    if which_sec_str == 3 | which_sec_str == 4
        turn = randi([0,1],1,1); %add or not beta-turn
        if turn == 1
	    for j=1:(fix(l(i)/2)-2)
		[phi, psi, name] = secStructs(which_sec_str);
                disp(name);
                potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
	    end
	    which_beta=randi([1,4],1,1);
            [phi, psi, name] = secStructs(2*which_beta+6);
            disp(name);
            potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
            [phi, psi, name] = secStructs(2*which_beta+7);
            disp(name);
            potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
	    for j=(fix(l(i)/2)+1):l(i)
		[phi, psi, name] = secStructs(which_sec_str);
                disp(name);
                potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
	    end
	else
	    for j=1:l(i)
		[phi, psi, name] = secStructs(which_sec_str);
                disp(name);
                potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
	    end
	end
    else
	for j=1:l(i)
            [phi, psi, name] = secStructs(which_sec_str);
            disp(name);
            potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
        end
    end
end

end
