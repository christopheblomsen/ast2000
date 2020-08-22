import numpy as np
import matplotlib.pyplot as plt


class integral:
    def __init__(self,Sigma,miu,a,b):
       self.Sigma=Sigma
       self.miu=miu
       self.a=a
       self.b=b

    def f(self,x):
        self.x=x 
        return 1 / (np.sqrt(2 * np.pi)*self.Sigma) * np.exp(-1/2*((x - self.miu) / self.Sigma)**2)

    def P(self):
        #integrasjon vha midtpunktsmetoden
        N=1000      #Antall tidssteg
        dx=(self.b-self.a)/N  #tidsintevall
        int=0
        x=self.a         #Start
        for i in range(N):

            int+=(self.f(x)+self.f((x+dx)))/2*dx
            x+=dx
        return int
