load tetrapepnew.pdb
select bb, name n+ca+c+o

hide all
show lines
show sticks, bb
center

set dash_gap, 0

distance /tetrapepnew//A/ALA`3/O, /tetrapepnew//A/LYS`4/CA
distance /tetrapepnew//A/LYS`4/O, /tetrapepnew//A/ASP`5/CA
distance /tetrapepnew//A/ASP`5/O, /tetrapepnew//A/VAL`6/CA

distance /tetrapepnew///UNK`9000/X1, /tetrapepnew///UNK`9000/X2
distance /tetrapepnew///UNK`9001/X1, /tetrapepnew///UNK`9001/X2
distance /tetrapepnew///UNK`9002/X1, /tetrapepnew///UNK`9002/X2

hide labels
