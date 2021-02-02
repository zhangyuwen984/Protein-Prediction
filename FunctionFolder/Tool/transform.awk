BEGIN {flag=0;FS="[,]" };
/xyz/ {flag=1; ic=0};
{   if (flag == 1)
    {   if ( $0 ~ /x/ && $0 ~ /y/ && $0 ~ /z/ && $0 !~ /xyz/  ) { print $1, $2, $3};
    }
}
