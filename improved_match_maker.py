## File name: match_maker.py
## Author: Alex Machtay (machtay.1@osu.edu)
## Revised by Dylan Wells (wells.1629@osu.edu) Summer 2022
## Date: 11/1/21
## Purpose:
##	This script is designed to read in the results from an XFdtd simulation. It will find the
##	 impedance information for the antenna designed by the GENETIS loop. The data is located
##	 inside the Run0001/output/SteadyStateOutput directory inside of a simulation number in 
##	 the XFdtd project directory.
#
# 
#
## Instructions:
##      To run, give the following arguments and a dog a bone
##      <source directory> <destination> <generation> <NPOP> <frequency>
## Example:
##      python match_maker.py ../properly_matched_antennas.xf . 0 50 14 
#
#
#
## Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import csv
import matplotlib.cm as cm
from scipy.optimize import minimize
import math
####import skrf as rf
####rf.stylely()
## Arguments
'''
parser = argparse.ArgumentParser()
parser.add_argument("source", help="Path to the Xfdtd project directory.", type=str)
parser.add_argument("destination", help="Location to output .cir files to", type=str)
parser.add_argument("gen", help="Generation of the run.", type=int)
parser.add_argument("NPOP", help="Population size.", type=int)
parser.add_argument("freq_num", help="Frequency number to match to (0-59)", type=int)
#parser.add_argument("num_files", help="Number of files to read through.", type=int)
g=parser.parse_args()

#
#
#
## Read in the files to obtain the impedances
file_base = "/Run0001/output/SteadyStateOutput/f"
impedances_R = [] # Real component of the impedance
impedances_I = [] # Imaginary component of the impedance
for i in range(g.gen*g.NPOP + 1, (g.gen+1)*g.NPOP + 1):
	## The numbering scheme of an XFdtd project always contains 6 digits, meaning it will have
	##  leading zeroes. We need to account for this so we read the correct files.
	if i < 10:
		lead = "00000" + str(i)
	elif i >=10 and i < 100:
		lead = "0000" + str(i)
	elif i >= 100 and i < 10**3:
		lead = "000" + str(i)
	elif i >= 10**3 and i < 10**4:
		lead = "00" + str(i)
	elif i >= 10**4 and i < 10**5:
		lead = "0" + str(i)
	else:
		lead = str(i)
	## Now we can read in the files and put the impedances into the lists defined above	
	with open(g.source + "/Simulations/" + lead + file_base + str(g.freq_num) + "/system.ssout", "r") as f:
		system_read = csv.reader(f, delimiter = " ")
		for line in f:
			if "Impedance" in line:
				impedances_R.append(float(line.split()[1]))
				impedances_I.append(float(line.split()[2][1:]))
				break ## need to break out of this because the phrase appears again in the next line
	f.close()
#
#

## Now we need to take the antenna impedances and the frequency to find the inductance and cpacitance
##  of the matching circuits.
## We'll assume the source resistance is 50 ohms
Rs = 50
freq = 10**6
#freq = (g.freq_num*16.67+83.33)*10**6
inductances = []
capacitances = []
inductancesHP = []
capacitancesHP = []
inductancesLP = []
capacitancesLP = []
## There are 8 different kinds of matching circuits that can be appropriate!!!
## See this explanation on L-network analysis to see why:
##  http://www.ittc.ku.edu/~jstiles/723/handouts/section_5_1_Matching_with_Lumped_Elements_package.pdf
## These 8 depend on the sign of the reactance and the magnitude of each component of 
##  the load compared to the source resistor
## See this, trying different values and signs: 
## https://www.analog.com/en/design-center/interactive-design-tools/rf-impedance-matching-calculator.html
##  (since the characteristic impedance is real, you can't get 2 of the designs)
#
#
#
## The 8 topologies are (first element is relative to the load/antenna):
##  series capacitor, parallel inductor
##  series inductor, parallel capacitor
##  series inductor, parallel inductor
##  series capacitor, parallel capacitor
##  parallel capacitor, series inductor
##  parallel inductor, series capacitor
##  parallel inductor, series inductor
##  parallel capacitor, series capacitor
#
#
#
## In practice, we'll only use 4 of these--the ones which don't repeat elements
## We will choose the first element to match the resistances and the second to match the impedances
## Here are the rules
##  If R < Rs (resistance of load less than source resistance):
##   Add element in series, add opposite element in parallel
##   To do this, we must match the *conductance* followed by *susceptance*
##  If R > Rs
##   Add element in parallel, add opposite element in series
##   To do this, we will match the *resistance* followed by *reactance*
#

'''
#
#
## Let's define functions to do this
## Naming Conventions are defined from the source to the load
## R < Rs; series capacitor, parallel inductor
def SLPC(Zreal, Zimag, Rs, f): # series inductor, parallel capacitor given impedance and characteristic resistance
	## We want to put an inductor in series such that the conductance matches that of Rs
	## Since we have a purely reactive element in series, the reistance of the series can't change
	## We have to match the conductance--this requires Re(Y)=1/Rs
	## Given these two conditions, we obtain:
	##  X = sqrt(Rs*Ra - Ra^2) - Xa (a for the antenna)
	## This tells us the reactance. We solve for the inductance (making sure it's positive!)
	Xl = np.sqrt(Zreal*Rs - Zreal*Zreal) - Zimag
	L = Xl/(2*np.pi*f)
	## Next, we need the capacitor
	## The reactance of the parallel element is found by matching the resistance
	Xc = Rs*complex(Zreal, Zimag + Xl)/(complex(Zreal - Rs, Zimag + Xl))
	C = abs(1/(2*np.pi*f*np.sqrt(-1 + 0j)*Xc))
	return C, L
#
#
#
## Let's repeat this for the other options
## Here, we flip C and L above
def SCPL(Zreal, Zimag, Rs, f): # series capacitor, parallel inductor given impedance and characteristic resistance
	## We want to put a capacitor in series such that the conductance matches that of Rs
	## Since we have a purely reactive element in series, the reistance of the series can't change
	## We have to match the conductance--this requires Re(Y)=1/Rs
	## Given these two conditions, we obtain:
	##  X = sqrt(Rs*Ra - Ra^2) + Xa (inductor)
	## This tells us the reactance. We solve for the inductance (making sure it's positive!)
	Xc = np.sqrt(Zreal*Rs - Zreal*Zreal) + Zimag ## note that the capacitor *adds* Z.imag
	
	C = 1/(2*np.pi*f*Xc)
	## Next, we need the inductor
	## The reactance of the parallel element is found by matching the resistance
	Xl = Rs*complex(Zreal, Zimag - Xc)/(complex(Zreal - Rs, Zimag - Xc)) ## note we subtract Xc
	L = abs(Xl/(2*np.pi*f*np.sqrt(-1 + 0j)))
	return C, L
#
#
#
## If R > Rs, we need to put the parallel component first
## Parallel inductor, series capacitor
def PLSC(Zreal, Zimag, Rs, f):
	## See the appendix in the paper!
	Xl = np.sqrt((Rs*Zreal**2 + Zimag**2*Rs-Rs**2*Zreal)/Zreal)
	L = (np.sqrt(Zreal*(Rs*Zreal**2 - Rs**2*Zreal + Zimag**2*Rs)) + Zimag*Rs)/(2*np.pi*f*(Zreal - Rs))
	Xc = -Xl
	C = abs(1/(2*np.pi*f*Xc))
	return C, L

#Parralel Capacitor, Series inductor
def PCSL(Zreal, Zimag, Rs, f):
        XC = (Zimag + np.sqrt(Zreal/Rs)*np.sqrt(Zreal**2+Zimag**2-Zreal*Rs))/(Zreal**2+Zimag**2)
        XL = ((1/XC) + (Zimag*Rs/Zreal) - Rs/(XC*Zreal))
        C = XC/(2*np.pi*f)
        L = XL/(2*np.pi*f)
        return C, L
	
#To broadband match we will cascade these L-matching elements
#NLC is the number of L matching circuits
#Need to create program that calculates NLC or define a standard
#Need to fix the functions to take in inputs for parts of impedance seperately, no Z class


#Creating a highpass and lowpass broadband matching

#Calculates the number of L circuits needed to match two resistances across a frequency range
def calcN(freqRange, Rl, Rs):
	BW = np.log(freqRange[1]/freqRange[0])/np.log(2)
	Q = ((2**BW)**(1/2))/(2**BW-1)
	N = abs(np.log(Rl/Rs)/np.log(1+Q**2))
	return int(math.ceil(N))


numCirc = calcN([100, 500], 10, 50) 
#Returns a list of capaitances and inductances to broadband match the circuit
def broadbandMatchHP(Rl, Xl, Rs, numCirc, freq):
	if Rs < Rl:
		#calculating the Q value from the number of L circuits we want
		Q = math.sqrt(math.e**(np.log(Rl/Rs)/numCirc)-1)	
		#we use this Q value to find the ratio between each imaginary resistance
		ratio = Q**2 + 1
		#creating a list of the real and imaginary resistances
		tempResistances = [Rs]
		for i in range(numCirc-1):
			tempResistances.append(tempResistances[-1]*ratio)
		tempResistances.append(Rl)
		tempCaps = []
		tempInducs = []
		#Calculating the capacitances and inductances to match each of these imaginary resistances
		for i in range(len(tempResistances)-1):		
			if tempResistances[i+1] == Rl:
				C, L = PLSC(Rl, Xl, tempResistances[i], freq)
				tempCaps.append(C)
				tempInducs.append(L)
			else:
				C, L = PLSC(tempResistances[i+1], 0, tempResistances[i], freq)
				tempCaps.append(C)
				tempInducs.append(L)
		return tempCaps, tempInducs
	else:
			
		Q = math.sqrt(math.e**(np.log(Rs/Rl)/numCirc)-1)	
		ratio = Q**2 + 1
		tempResistances = [Rl]
		for i in range(numCirc-1):
			tempResistances.append(tempResistances[-1]*ratio)
		tempResistances.append(Rs)
		tempCaps = []
		tempInducs = []
		for i in range(len(tempResistances)-1):
			if i == 0:
				C, L = SCPL(Rl, Xl, tempResistances[i+1], freq)
				tempCaps.append(C)
				tempInducs.append(L)
			else:
				C, L = SCPL(tempResistances[i], 0, tempResistances[i+1], freq)
				tempCaps.append(C)
				tempInducs.append(L)
		tempCaps.reverse()
		tempInducs.reverse()
		return tempCaps, tempInducs
#lowpass version of the broadbandmatch, either version can be used
def broadbandMatchLP(Rl, Xl, Rs, numCirc, freq):
        if Rs < Rl:
                #calculating the Q value from the number of L circuits we want
                Q = math.sqrt(math.e**(np.log(Rl/Rs)/numCirc)-1)
                #we use this Q value to find the ratio between each imaginary resistance
                ratio = Q**2 + 1
                #creating a list of the real and imaginary resistances
                tempResistances = [Rs]
                for i in range(numCirc-1):
                        tempResistances.append(tempResistances[-1]*ratio)
                tempResistances.append(Rl)
 #               print(tempResistances, numCirc)
                tempCaps = []
                tempInducs = []
                #Calculating the capacitances and inductances to match each of these imaginary resistances
                for i in range(len(tempResistances)-1):
                        if tempResistances[i+1] == Rl:
                                C, L = PCSL(Rl, Xl, tempResistances[i], freq)
                                tempCaps.append(C)
                                tempInducs.append(L)
                        else:
                                C, L = PCSL(tempResistances[i+1], 0, tempResistances[i], freq)
                                tempCaps.append(C)
                                tempInducs.append(L)
                return tempCaps, tempInducs
        else:

                Q = math.sqrt(math.e**(np.log(Rs/Rl)/numCirc)-1)
                ratio = Q**2 + 1
                tempResistances = [Rl]
                for i in range(numCirc-1):
                        tempResistances.append(tempResistances[-1]*ratio)
                tempResistances.append(Rs)
#                print(tempResistances, numCirc)
                tempCaps = []
                tempInducs = []
                for i in range(len(tempResistances)-1):
                        if i == 0:
                                C, L = SLPC(Rl, Xl, tempResistances[i+1], freq)
                                tempCaps.append(C)
                                tempInducs.append(L)
                        else:
                                C, L = SLPC(tempResistances[i], 0, tempResistances[i+1], freq)
                                tempCaps.append(C)
                                tempInducs.append(L)
                return tempCaps.reverse(), tempInducs.reverse()
#Transforming the broadband LP match to a broadband bandpass

def bandPass(Rl, Xl, Rs, numCirc, freq):
	C, L = broadbandMatchLP(Rl, Xl, Rs, numCirc, freq)
	if Rs < Rl:
		Q = math.sqrt(math.e**(np.log(Rl/Rs)/numCirc)-1)
	else:
		Q = math.sqrt(math.e**(np.log(Rs/Rl)/numCirc)-1)
	Parallels = []
	Series = []
	for ele in C:
		Parallels.append([Q*ele, 1/(Q*ele)])
	for ele in L:
		Series.append([Q*ele, 1/(Q*ele)])
	return Parallels, Series

N = calcN([100, 1000], 229.839, 50)
C, L = broadbandMatchLP(229.839, -151.515, 50, 14, 316.57)
print(N,"\n", C)
print("\n")
print(L)

'''
freqRange = [100*10**6, 1000*10**6]
for i in range(len(impedances_R)):
	N = calcN(freqRange, impedances_R[i], Rs)
	C, L = broadbandMatchHP(impedances_R[i], impedances_I[i], Rs, N, np.sqrt(freqRange[0]*freqRange[1]))
	#C, L = broadbandMatchHP(311, -197, Rs, 1, 300*10**6)
	capacitances.append(C)
	inductances.append(L)	
#print("Antenna impedance and frequency (MHz)")
#print(impedances_R[48], impedances_I[48],"j", freq/10**6)
#print("Calculated inductance (nH) and capacitance (pF)")
#print(inductances[48]/10**-9, capacitances[48]/10**-12)
## Right now, we are using a three port system, with an inductor in parallel and a capacitor
##  in series.

for i in range(g.gen*g.NPOP + 1, (g.gen+1)*g.NPOP + 1):
	with open(g.destination + "/" + str(i) + "_matching_circuit.cir", "w") as f:
		f.writelines(["TEST","\n"])
		f.writelines("\n")
		teststring = ""
		for j in range(len(capacitances[i-1])+2): #Generating a string of nodes needed
			teststring = teststring + " P" + str(j+1)
		f.writelines([".SUBCKT match_test "+teststring,"\n"]) #.SUBCKT starts the netlist for a subcircuit in Xf
		if impedances_R[i-1] > Rs: #Downward Highpass Match 
			for j in range(len(capacitances[i-1])):
				print(len(capacitances[i-1]))
				if j == 0: #Xfdtd had the 'parallel' components connect to P2 in the example
					f.writelines(["C01 P1 P3 ", "{:e}".format(capacitances[i-1][j]),"\n"])
					f.writelines(["L01 P3 P2 ", "{:e}".format(inductances[i-1][j]),"\n"])
				else:
					f.writelines(["C0"+str(j+1)+" P"+str(j+2)+" P"+str(j+3)+" ", "{:e}".format(capacitances[i-1][j]),"\n"])
					f.writelines(["L0"+str(j+1)+" P"+str(j+3)+" P2 ", "{:e}".format(inductances[i-1][j]),"\n"])
		else: #Upward Highpass Match
			for j in range(len(capacitances[i-1])):
				print(len(capacitances[i-1]))
				if j == 0:
					f.writelines(["L01 P1 P2 ", "{:e}".format(inductances[i-1][j]),"\n"])
					f.writelines(["C01 P1 P3 ", "{:e}".format(capacitances[i-1][j]),"\n"])
				else:
					f.writelines(["L0"+str(j+1)+" P"+str(j+2)+" P2 ", "{:e}".format(inductances[i-1][j]),"\n"])
					f.writelines(["C0"+str(j+1)+" P"+str(j+2)+" P"+str(j+3)+" ", "{:e}".format(capacitances[i-1][j]),"\n"])
		f.writelines([".ENDS match_test","\n"])
		f.writelines("\n")
		f.writelines(".END")
	f.close()

'''

