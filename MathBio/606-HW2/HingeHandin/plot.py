#! /usr/bin/env python
import sys
import pylab as plt

def readfile(str):
	#Initialize three arrays to capture values of angle, Energy and nltcut
	a=[]
	b=[]
	c=[]
	
	fh=open("interdomain_interaction_hinge_models_res"+str+".txt", 'r')
	for line in fh:
		if line[0]=='#':
			continue

		(ang,cal,ngt,nlt)=line.split(' ')
		#print line.split(' ')
		a.append(float(ang))
		b.append(float(cal))
		c.append(float(nlt))
	#print angles
	#print coul
	#print cutoff
	return (a,b,c)

#Populate data into these data_resX arrays so that they can be plotted later

if len(sys.argv)!=4:
	sys.exit("Usage is python plot.py res1 res2 res3")

(angle_res1, energy_res1, nlt_res1)= readfile(sys.argv[1])
# HOMEWORK: READ IN FILES FOR RESIDUES 2 AND 3. Hint: check the next lines to find variable names to use

(angle_res1, energy_res1, nlt_res1) = readfile(sys.argv[1])
(angle_res2, energy_res2, nlt_res2) = readfile(sys.argv[2])
(angle_res3, energy_res3, nlt_res3) = readfile(sys.argv[3])


#Identify maximum energy and nlt in order to figure out cutoffs for plotlabels
maxcaE= max( max(energy_res1), max(energy_res2), max(energy_res3))
maxnlt= max( max(nlt_res1), max(nlt_res2), max(nlt_res3))


#Plot for Calcium EE
plot1= plt.subplot(211) # Make 2 subplots, assign this plot to first column, first row

plot1.scatter( angle_res1, energy_res1, s=40, c='r', marker='o', edgecolors='none')
plot1.plot(angle_res1, energy_res1,c='r')

plot1.scatter(angle_res2, energy_res2, s=40, c='b', marker='o', edgecolor='none' )
plot1.plot(angle_res2, energy_res2, c='b')

plot1.scatter(angle_res3, energy_res3, s=40, c='y',marker='o', edgecolor='none' )
plot1.plot(angle_res3, energy_res3, c='y')

plot1.set_ylim(90, maxcaE+20) # what we found before
plot1.set_xticks([0,60,120,180,240,300,360]) #Where to place markers
plt.xlabel("phi(degrees)")
plt.ylabel('Calcium EE, Kcal/mol')
plot1.legend( [sys.argv[1],sys.argv[2],sys.argv[3]], loc=4)

#HOMEWORK: COMPLETE PLOTS FOR NLT
plot2=plt.subplot(212) #2 plots, first column second row

plot2.scatter( angle_res1, nlt_res1, s=40, c='r', marker='o', edgecolors='none')
plot2.plot(angle_res1, nlt_res1,c='r')

plot2.scatter(angle_res2, nlt_res2, s=40, c='b', marker='o', edgecolor='none' )
plot2.plot(angle_res2, nlt_res2, c='b')

plot2.scatter(angle_res3, nlt_res3, s=40, c='y',marker='o', edgecolor='none' )
plot2.plot(angle_res3, nlt_res3, c='y')

plot2.set_ylim(0,maxnlt+20)
plot2.set_xticks([0,60,120,180,240,300,360])
plt.xlabel("phi(degrees)")
plt.ylabel('Number of clashes')
plot2.legend( [sys.argv[1],sys.argv[2],sys.argv[3]], loc=4)


#plt.show()
#This command saves the graph generated
plt.savefig('homeworkplot2.png')
#show()
