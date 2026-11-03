#! /usr/bin/env python

# import std librairies

import numpy as np
import pylab as plt

def smoothmin_bad(a,b,p):
    return 0.5*(a+b-smoothabs_bad(a-b,p))

def smoothmin_good(a,b,p):
    return 0.5*(a+b-smoothabs_good(a-b,p))

def smoothmin_sharp(a,b):
    return 0.5*(a+b-abs(a-b))

def smoothabs_bad(c,p):
    sf=1.0/p
    rsf=1.0/sf
    absc=np.zeros((len(c)))
    for i in range(len(c)):
        if c[i]<-p:
            absc[i]=-c[i]
        elif c[i]>p:
            absc[i]=c[i]
        else:
            absc[i]=np.log(0.5*(np.exp(c[i]*sf)+np.exp(-c[i]*sf)))*rsf
    return absc

def smoothabs_good(c,p):
    if p!=0.0:
        sf=1.0/p
    else:
        sf=1.0
    rsf=1.0/sf
    absc=np.zeros(np.shape(c))
    absc=np.log(np.exp(c*sf)+np.exp(-c*sf))*rsf

    return absc

if __name__ == "__main__":

    x=np.linspace(-20,20,1000)

    # sx_bad=smoothmin_bad(0.0,x,20)
    sx_good_100=smoothmin_good(0.0,x,100)
    sx_good_20=smoothmin_good(0.0,x,10)
    sx_sharp=smoothmin_sharp(0.0,x)

    plt.figure(1)
    plt.plot(x,sx_good_100,'b-', label='100')
    plt.plot(x,sx_good_20,'r-', label='20')
    plt.plot(x,sx_sharp,'g-', label='0')
    plt.grid()
    plt.legend()
    plt.show()