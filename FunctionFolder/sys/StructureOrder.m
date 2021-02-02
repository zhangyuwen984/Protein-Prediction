function S_order = StructureOrder(f, volume, numIons, deltaFing, weight)

[r, c] = size(f);
order_full = 0;
for j = 1 : r
 orderAB = 0;
 for k = 1 : c
   orderAB = orderAB + f(j,k)^2;
 end
 orderAB = orderAB*deltaFing/power(volume/sum(numIons),1/3);
 order_full = order_full + orderAB*weight(j);
end
S_order = sqrt(order_full);
