BEGIN {flag=0}; 
/Final stress tensor components/ {flag=1; ic=0}; 
{if (flag == 1) 
    {ic++; 
	if (ic == 4) {p11=$2; p23=$4; p32=$4}; 
	if (ic == 5) {p22=$2; p13=$4; p31=$4}; 
	if (ic == 6) {p33=$2; p12=$4; p21=$4}; 
    }
};
END {print p11*(-10),p12*(-10),p13*(-10); print p21*(-10),p22*(-10),p23*(-10); print p31*(-10),p32*(-10),p33*(-10)}
