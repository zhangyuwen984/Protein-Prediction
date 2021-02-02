/in kB/ {
    p11=$3; p22=$4; p33=$5; p12=$6; p23=$7; p13=$8;
    { print p11, p12, p13;
      print p12, p22, p23;
      print p13, p23, p33 }
}
