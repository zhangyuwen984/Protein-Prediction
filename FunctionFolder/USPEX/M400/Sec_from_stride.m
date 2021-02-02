function [names, angles_list] = Sec_from_stride(pdb, sec_len)

angles_list=[];
all_angles_list=[];
%[stat, stride_out] = unix([ORG_STRUC.USPEXPath '/FunctionFolder/USPEX/M400/stride ' pdb]);
[stat, stride_out] = unix(['stride ' pdb]);
    k = strfind(stride_out, 'ASG ');
    if stat == 0
        for i=1:(size(k, 2))
            if i < (size(k, 2))
                sec_string = stride_out(k(i):k(i+1));
            else
                sec_string = stride_out(k(end):end);
            end
            sec_structure = textscan(sec_string, '%s %s %s %d %d %s %s %f %f %f %s');
            all_angles_list = [ all_angles_list; [sec_structure{8} sec_structure{9}]];
	    try
        	all_names{i}=sec_structure{7}{1};
	    catch ME
         	causeException = MException('MATLAB:file', pdb);
         	ME = addCause(ME,causeException);
               rethrow(ME)
	    end
        end
	size(all_angles_list, 1);
	sec_len;
	k = randi([1, size(all_angles_list, 1)-sec_len+1],1,1);
        angles_list=all_angles_list(k:k+sec_len-1,:);
	for j=1:sec_len
           names{j}=all_names{k+j-1};
	end
    else
        which_sec_str=randi([1,7],1,1);
	names=[];
        for i=1:sec_len
	    [phi, psi, name] = secStructs(which_sec_str);
            names = [names; name];
            angles_list = [angles_list; [double(phi), double(psi)]];
        end
    end
end
