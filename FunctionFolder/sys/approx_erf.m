% determines erf(datapoints) using erf function table : from -4.01 to 4.01 with step 0.01
function [my_erf]=approx_erf(datapoints, erf_table)
len=length(datapoints);
my_erf = zeros(len,1);
for i=1:len
    if datapoints(i)>=4
        my_erf(i)=1;
    elseif datapoints(i)<=-4
        my_erf(i)=-1;
    else
        x1 = floor(datapoints(i)/0.01+402);
        y3 = erf_table(x1+2);
        y2 = erf_table(x1+1);
        y1 = erf_table(x1);
        y0 = erf_table(x1-1);
        a0 = -y0/6 + y1/2 - y2/2 + y3/6;
        a1 = y0/2 - y1 + y2/2;
        a2 = -y0/3 - y1/2 + y2 - y3/6;
        a3 = y1;
        q = (datapoints(i)/0.01+402-x1);
        my_erf(i) = ((a0*q+a1)*q+a2)*q+a3;
    end
end
