# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 21:59:46 2014

@author: sophie
"""

x = raw_input ("Input x:  ")
y = raw_input ("Input y:  ")


x = int(float(x))
y = int(float(y))
   
if x > y:
    print '1'
elif x == y:
    print 'No, not right'
elif x < y:
    print '-1'