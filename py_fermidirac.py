# -*- coding: utf-8 -*-
import numpy as np

def fermidirac(x,s):

    f = lambda k: np.sqrt(x**2+np.pi**2*(2*k-1)**2)

    if s==0.5:    # F_{1/2}(x)
        a = np.array((1.0/770751818298,-1.0/3574503105,-13.0/184757992,85.0/3603084,3923.0/220484,74141.0/8289,-5990294.0/7995))
        g = lambda k:np.sqrt(f(k)-x)
    
    else:         # F_{-1/2}(x)

        a = np.array((-1.0/128458636383,-1.0/714900621,-1.0/3553038,27.0/381503,3923.0/110242,8220.0/919))
        g = lambda k:-0.5*np.sqrt(f(k)-x)/f(k)

    F = np.polyval(a,x) + 2*np.sqrt(2*np.pi)*sum(map(g,range(1,21)))
    return F


if __name__ == '__main__':

    n = 11
    x = np.linspace(-6,10,n)

    f = np.frompyfunc(lambda s:fermihalf(s,1),1,1)
    print(f(x))












