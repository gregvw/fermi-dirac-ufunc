import numpy as np
import fermidirac as fd 
import py_fermidirac as pyfd
import timeit as ti

if __name__ == '__main__':

    m = 101
    x = np.linspace(-6,10,m)
    f = fd.half(x)
    print(f)  
   
