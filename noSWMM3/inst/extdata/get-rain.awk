/TIMESERIES/ { partof=1 }
/REPORT/ { partof=0 }
{ if ( (partof==1) && ($2 ~ /FILE/) )  {print $3} }