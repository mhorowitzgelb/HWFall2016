$title HW1-1
$offsymxref offsymlist offuelxref offuellist offupper
option limrow=0, limcol=0;

free variable obj "total loss";

positive variables
x1,
x2,
x3;

equations
loss,
eq1,
eq2,
eq3;

loss..
obj =e= 3*x1 + 2*x2 - 33*x3;

eq1..
x1 - 4*x2 + x3 =l= 15;

eq2..
9*x1  + 6*x3 =e= 12;

eq3..
-5 * x1 + 9 * x2 =g= 3;

model prob1 /all/;

solve prob1 using lp minimize obj;

parameter x1val, x2val, x3val, objval;

objval = obj.l;
x1val = x1.l;
x2val  = x2.l;
x3val = x3.l;

display objval, x1val, x2val, x3val;



