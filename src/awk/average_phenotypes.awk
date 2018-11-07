BEGIN {printf "accession_id\tphenotype_value\n"} 
{a[$1] = a[$1]+$2; b[$1] = b[$1]+1} 
END {for(i in a) print i "\t" a[i]/b[i]}
