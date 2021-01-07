#-----------------------------------------------------------------------------
#PHY3042 - Assignment 1
#Student: 6463360
#To be viewed with the PDF for reference to equations.
#Written in Spyder (Python 3.8)
#-----------------------------------------------------------------------------


#---Import all useful modules here--------------------------------------------
import math as m
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


#---Variables that will be used are defined here------------------------------
#Current values are recommended.
x = 5               #The range of x-values = -x..x.
dx = 0.1            #The change from x..x1.
xrange = int(x/dx)  #The number of x-values in a set.
t = 10              #The range of t-values = 0..t.
dt = 0.1            #The change from t..t1.
trange = int(t/dt)  #The number of t-values in a set.

psi = np.empty(2*xrange, complex)            #This will be 'current' psi.
psiplus = np.empty(2*xrange, complex)        #This will be 'next' psi.
exppsi = np.empty(trange, complex)           #The expectation value of x.
expT = np.empty(trange)                      #Holds values of t to trange.
allT = np.empty(2*xrange)                    #This is t-values increasing.
allX = np.empty(2*xrange)                    #This is x-values increasing.
A = np.empty(2*xrange, complex)              #A array from eq. 6.7, calculated
                                             #using eq. 6.4.
B = np.empty(2*xrange, complex)              #B array from eq. 6.7 calculated
                                             #using eq. 6.5.
U = np.empty(2*xrange, complex)              #U array from eq. 7.11 and 7.12.
R = np.empty(2*xrange, complex)              #R array from eq. 7.21 and 7.22.
#A temp variable will appear throughout. This will be initialised on each
#use so it best suits it's requirement.


#---Part a--------------------------------------------------------------------
#---Calculation of Normalised psi---------------------------------------------
def nopsi(k):                                   #Calculates psi.
    return m.exp(-(k-1)**2)
def nopsi_conj(k):                              #Calculates psi-conjugate.
    return np.conjugate(m.exp(-(k-1)**2))
def nopsi_sq(k):                                #Calculates psi-squared.
    return nopsi(k)*nopsi_conj(k)


N, err = quad(nopsi_sq, float('-inf'), float('inf'))
N = m.sqrt(1/N)


def Npsi(k, l):                                  #Calculates psi with N.
    return l*nopsi(k)


#---Make t=0 array------------------------------------------------------------
for j in range(-xrange, xrange):        #Populates psiplsu with t=0 values.
    psiplus[j+xrange] = Npsi(j*dx, N)


#---Make the A array from eq. 6.4---------------------------------------------
#Due to the nature of the equations that will be done in eq. 7.1 and eq. 7.2
#The A array can be 1 dimentional and hold only the Aj values.
#This makes the following much more simple than if a multidimentional array
#were to be used.
for j in range(-xrange, xrange):
    A[j+xrange] = complex(-2.0 - (dx**2)*((j*dx)**2), (4.0*(dx**2))/dt)


#---Initialise ALL array------------------------------------------------------
psi = psiplus


#---Populate allT array with values-------------------------------------------
temp1 = np.ones(2*xrange)
allT = 0*temp1
for tn in range(1, trange):
    allT = np.concatenate((allT, tn*temp1))


#---Populate allX array with values-------------------------------------------
temp2 = np.empty(2*xrange)
for j in range(-xrange, xrange):
    temp2[j+xrange] = j*dx
for j in range(trange):
    allX[j] = temp2[j]
for j in range(trange-1):
    allX = np.concatenate((allX, temp2))


#---Calculate expectation value at t=0----------------------------------------
#This is the trapesium rule method of calculating <x>.
#temp3 = np.empty(2*xrange, complex)
#---Calculate exp-value of x at tn---
#for xn in range(-xrange, xrange):    #Function from expectation equation.
#    j = xn + xrange
#    temp3[j] = xn*dx * np.conjugate(psiplus[j]) * psiplus[j]

#s = 0                                #Initialise sum at 0.
#for j in range(1, 2*xrange-2):       #Sum of temp[1] to temp[2*xrange-1].
#    s = s + psiplus[j]

#Use trapesium rule to calculate expectation value for x at time tn.
#exppsi[0] = 0.5*dx*(temp3[0] + temp3[2*xrange-1] + (2*s))
#expT[0] = 0

#This is the averaging method of calculating <x>.
exppsi[0] = np.average(psiplus)
expT[0] = 0

#---Loop over trange to calculate psi plus at all t-values--------------------
for tn in range(trange-1):


    #---Create B using eq. 6.5---
    for xn in range(-xrange+1, xrange-1):
        i = xn + xrange
        B[i] = - psiplus[i+1] - psiplus[i-1] + \
            psiplus[i]*complex(2.0 + (dx**2)*((xn*dx)**2), 4.0*(dx**2)/dt)


    #---Make U and R required for eq .7.4---
    U[0] = 1/A[0]                   #First U-value.
    for j in range(2*xrange-1):
        U[j+1] = 1/(A[j+1] - U[j])

    R[0] = B[0]*U[0]                #First R-value.
    for j in range(2*xrange-1):
        R[j+1] = (B[j+1] - R[j])*U[j+1]


    #---Calculate the psiplus array---
    psiplus[2*xrange-1] = R[2*xrange-1]     #Final value of psiplus.
    for j in reversed(range(2*xrange-2)):
        psiplus[j] = R[j] - U[j]*psiplus[j+1]


    #---Make and concatenate psi array---
    psi = np.concatenate((psi, psiplus))


    #---Calculate exp-value of x at tn---
#This is the trapesium rule method of calculating <x>.
#    for xn in range(-xrange, xrange):    #Function from expectation equation.
#        j = xn + xrange
#        temp3[j] = xn*dx * np.conjugate(psiplus[j]) * psiplus[j]

#    s = 0                                #Initialise sum at 0.
#    for j in range(1, 2*xrange-2):       #Sum of temp[1] to temp[2*xrange-1].
#        s = s + psiplus[j]

    #Use trapesium rule to calculate expectation value for x at time tn.
#    exppsi[tn+1] = np.average(psiplus)
#    expT[tn+1] = tn*dt
    
    
    #This is the averaging method of calculating <x>.
    exppsi[tn+1] = np.average(psiplus)
    expT[tn+1] = tn*dt


#---Make 3 by 2*xrange*trange array of all values---
#Stacks the psi,  x and t values into one array for plotting.
allpsi = np.stack((allX, psi, allT), axis=1)
#Make expeectation value array plottable.
expect = np.stack((exppsi, expT), axis=1)


#---Now to plot the data------------------------------------------------------
#Plots the <x> as a function of t.
plt.title("<x> as a function of t")
plt.xlabel("t")
plt.ylabel("<x>")
plt.plot(expect[:,1].real, expect[:,0].real, color="green")
plt.show()


#Saves the psi(x,t) data to a text file.
#This is used to plot a 3D graph in GNUplot.
#A space between datasets of tn is required to plot the surface in GNUplot
#which is easily done with some trivial function in a text editor.
np.savetxt('psi.txt', allpsi.real)
