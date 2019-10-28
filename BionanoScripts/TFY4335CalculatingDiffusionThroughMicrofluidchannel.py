import pandas
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy
import numpy as np

l=1559.05 #micrometers
v=13/60 #microliters per second
h=500 #micrometers
w=75 #micrometers
t = l*h*w/v*10**(-9)
k_b=1.38*10**(-23)
T=293
visc=1.002*10**(-3)
pi=3.14
print(t)
plt.style.use("ggplot")

def func(z,a,b,c):
   return a*scipy.special.erf(z/(10**6*b)) + c

#Import .csv files, where the intensity values from ImageJ are stored
df1 = pandas.read_csv('FirstPointOutletGrayValue.csv')
df2 = pandas.read_csv('SecondPointOutletGrayValue.csv')

#create panda series from panda dataframe
firstXs= df1["X"]
secondXs = df2["X"]
firstYs= df1["Y"]
secondYs = df2["Y"]


#convert to list and make same size
firstXlist= firstXs.tolist()
secondXlistTemp= secondXs.tolist()
secondXlist = secondXlistTemp[:len(firstXlist)]
firstYlist= firstYs.tolist()
secondYlistTemp= secondYs.tolist()
secondYlist=secondYlistTemp[:len(firstXlist)]


#Adjust x-values for x=0 at center
x1Values = [-(i-firstXlist[-1]/2) for i in firstXlist]
x2Values = [-(i-secondXlist[-1]/2) for i in secondXlist]

# scale for consentration, since intensity increases as conserntration decreases
#Normalize, intensity for c0 is 105.67, intenisty for c=0 is 145.9
c1 = [((np.log10(145.90)-np.log10(i))/(np.log10(145.90)-np.log10(105.67))) for i in firstYlist]
c2 = [((np.log10(145.90)-np.log10(i))/(np.log10(145.90)-np.log10(105.67))) for i in secondYlist]

params, extras = curve_fit(func, x1Values, c1,p0=[0.5,10**-5,0.5])
params2, extras2 = curve_fit(func, x2Values, c2,p0=[0.5,10**-5,0.5])

plt.figure()
plt.plot(x1Values,c1, "k-", alpha=0.5)
plt.plot(x1Values,func(x1Values, *params), "r-", alpha=0.5)
plt.xlabel("Position [μm]")
plt.ylabel("Normed concentration")
print(params)
print(params2)
plt.savefig("Concentration_distribution_and_fitted_error_funnction_first_measurement.jpeg")

plt.figure()
plt.plot(x2Values,c2, "k-", alpha=0.5)
plt.plot(x2Values,func(x2Values, *params2), "r-", alpha=0.5)
plt.xlabel("Position [μm]")
plt.ylabel("Normed concentration")
plt.savefig("Concentration_distribution_and_fitted_error_funnction_second_measurement.jpeg")
def diffusionConstant(d1,d2,dt):
    return(d2**2-d1**2)/(4*dt)
    
D= diffusionConstant(params[1],params2[1],t)
print("The calculated diffusion-constant is: ")
print(D)
print("Which gives a particle size of ")

def size(D):
    return k_b*T/(6*pi*visc*D)*2
print(size(D))
print("The variance of the function is")
var1=np.sqrt(np.diag(extras)[1])
var2=np.sqrt(np.diag(extras2)[1])
def variance(d_1,d_2):
    return (2*(d_1**4)+2*(d_2**4))/(16*t**2)

print(variance(var1,var2))
print(np.sqrt(variance(var1,var2)))