function magString=magTypeString(mag)

switch mag
  case  1
    magString='  NM  ';
  case  2
    magString=' FM-HS';
  case -2
    magString=' FM-LS';
  case  3
    magString=' AFM-H';
  case -3
    magString=' AFM-L';
  case  4
    magString=' FM-LH';
  case  5
    magString=' AF-LH';
end

