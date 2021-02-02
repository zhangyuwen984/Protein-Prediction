function sec_parts = combineSecStructs(res_num)
% $Rev: 586 $
% $Author: maxim $
% $Date: 2014-08-30 11:58:30 +0400 (Sat, 30 Aug 2014) $

while 1
    a = RandInt(1,7,[0 res_num]);
    sec_parts = round(a/sum(a)*res_num);
    if sum(sec_parts) == res_num
        break
    end
end
