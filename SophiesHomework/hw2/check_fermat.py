# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 21:08:15 2014

@author: sophie
"""



#def check_fermat(a,b,c,n):

a = raw_input ("Input A:  ")
b = raw_input ("Input B:  ")
c = raw_input ("Input C:  ")
n = raw_input ("Input N greater than 2:  ")    

a = int(float(a))
b = int(float(b))
c = int(float(c))
n = int(float(n))

x=(a**n) + (b**n)
y=(c**n)
    
if x == y:
    print 'Holy smokes, Fermat was wrong!'
else:
        print 'No, not right'