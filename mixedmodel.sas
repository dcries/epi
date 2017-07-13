proc import datafile= "C:/Users/dcries/data/meas7.csv" 
	replace
	out=mydata;
  getnames=yes;
run;


proc mixed;
	class rep dow id;
	model modvigmin = rep dow;
	random id;
run;

proc mixed;
	class rep dow id;
	model modvigmin = rep dow;
	random id;
	repeated rep / subject=id type=ar(1);
run;
