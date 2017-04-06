NF<6 { barrier= $3/$2 ; printf "%-11s     1 %6.3f   %5.1f  %3.0f\n", $1, barrier, $4, $5}
NF>5 { barrier= $3/$2 ; pos = index($0, $6) ;
printf "%-11s     1 %6.3f   %5.1f  %3.0f %s\n", $1, barrier, $4, $5, substr($0,pos) }
