# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 20:54:07 2014

@author: sophie
"""

print 'Question 3.5'

def plusbar():
    print '+ - - - -',

def bar():
    print '|        ',

def twice(f):
    f()
    f()
    
def four(f):
    twice(f)
    twice(f)

def eight(f):
    four(f)
    four(f)

def plus_row():
    twice(plusbar)
    print '+'
    
def bar_row():
    twice(bar)
    print '|'
    
def gridhalf():
    plus_row()
    four(bar_row)

def grid():
    twice(gridhalf)
    plus_row()
    
grid()

def plus_row2():
    four(plusbar)
    print '+'
    
def bar_row2():
    four(bar)
    print '|'
    
def gridhalf2():
    plus_row2()
    eight(bar_row2)

def grid2():
    twice(gridhalf2)
    plus_row2()
    
grid2()