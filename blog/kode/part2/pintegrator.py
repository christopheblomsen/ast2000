# Egen kode
import numpy as np
import matplotlib.pyplot as plt
import math

class pintegrator():

    def __init__(self,f):
        self.f =f


    def integrate(self,r0 ,v0, T, dt):
        pass


class EulerCromer(pintegrator):

    def __init__(self,f):
        pintegrator.__init__(self,f)

    def integrate(self,r0,v0, T, dt):
        N = int(T//dt)
        r = np.zeros(N,float)
        t = np.linspace(0,T,N,dtype=float)
        r[0] = r0
        v = v0
        for i in range(N-1):
            ai = self.f(r[i],t[i])
            print(f'ai {ai}')
            v += ai*dt
            r[i+1] = r[i] + v*dt

        return r, t

class LeapFrog(pintegrator):

    def __init__(self,f):
        pintegrator.__init__(self,f)

    def integrate(self,r0,v0, T, dt):
        N = int(T//dt)

        if type(r0) is int:
            r = np.zeros(N,float)
        else:
            r = np.zeros((N,r0.shape[0],r0.shape[1]),float)
        t = np.linspace(0,T,N,dtype=float)

        r[0] = r0
        v = v0

        ai = self.f(r[0],t[0])
        for i in range(N-1):
            r[i+1] = r[i] + v*dt + 0.5*ai*dt**2
            ai_pluss1 = self.f(r[i+1],t[i+1])
            v += 0.5*(ai + ai_pluss1)*dt
            #print(f'v {v} | {np.sqrt(v[0]**2+v[1]**2)}')
            ai = ai_pluss1

        return r, t, v


if __name__ == '__main__':

    m = 1.5 #Vekt i kg
    def A1(r,t):
        return -9.8*m
    def A2(r,t):
        return -np.sin(r) # Beregn aksellerasjonen
    def A3(t):
        return t

    ecromer = EulerCromer(A2)
    #x2,t2 = ecromer.integrate(0,40,6,0.2)
    x2,t2 = ecromer.integrate(1,0,4*np.pi,0.01)
    plt.plot(t2,x2,label="EulerCromer")
    leapfrog = LeapFrog(A2)
    #x3,t3,v = leapfrog.integrate(0,40,6,0.2)
    x3,t3,v = leapfrog.integrate(1,0,4*np.pi,0.0001)
    #plt.plot(t3,x3,label="LeapFrog")

    t4 = np.arange(0,4*np.pi,0.001)
    x4 = np.cos(t4)
    plt.plot(t4,x4,label="Exact")
    plt.legend()
    plt.show()
