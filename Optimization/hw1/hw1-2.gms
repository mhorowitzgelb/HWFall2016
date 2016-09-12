$title HW1-1
$offsymxref offsymlist offuelxref offuellist offupper
option limrow=0, limcol=0;

SETS
	J /1,2,3/;


free variable obj;

positive variable x(J);

equations
	gain,
	eq1,
	eq2(J);

gain ..
	obj =e= 5 * (x("1") + 2 * x("2")) - 11 * (x("2") - x("3"));

eq1 ..
	3*x("1") =g= sum(J, x(J));

eq2(J) ..
	x(j) =l= 3;

model prob2 /all/;
solve prob2 using lp maximizing obj;	
	
display x.l, x.lo, x.up, prob2.objval;
