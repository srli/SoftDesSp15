# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:34:57 2014

@author: sophie li i am artist
collaboration with Zoher and Jay Woo

No unit tests with this assignment, since the things we generate are random.
"""
from random import randint
import Image
from math import *

def build_random_function(min_depth, max_depth):
    """Takes minimum function depth and maximum function depth as inputs, and
    generates a function using a set of 7 "blocks".
    """
    hello = ['x','y']
    func = ['x','y','cos_pi','sin_pi','prod','square','average']
    if max_depth == 1:
        return hello[randint(0,1)]
    else:
        block = func[randint(2,6)]
        if block == 'prod' or 'average': #accouts for when a block requires two inputs
            return [block, build_random_function(min_depth-1, max_depth-1), build_random_function(min_depth-1, max_depth-1)]
        elif not block == 'prod':
           return [block, build_random_function(min_depth-1, max_depth-1)]

def evaluate_random_function(f, x, y):
    """Takes function as input and a numerical value for x and y.
    This function evalues the input function using the values provided for x and y.
    Since f is length two unless the block to be evaluated requires two inputs,
    we can do calculations having our if statements look at the first statement in the list.
    """
    if f[0] == 'x': #If the first index is x or y, we've already reached the innermost layer and can stop our recursion
        return x
    elif f[0] == 'y':
        return y
    elif f[0] == 'square':
        return evaluate_random_function(f[1],x,y)**2
    elif f[0] == 'average':
        return (evaluate_random_function(f[1],x,y)+evaluate_random_function(f[2],x,y))/2
    elif f[0] == 'cos_pi':
        return cos(pi*evaluate_random_function(f[1],x,y))
    elif f[0] == 'sin_pi':
        return sin(pi*evaluate_random_function(f[1],x,y))
    elif f[0] == 'prod':
        return evaluate_random_function(f[1],x,y)*evaluate_random_function(f[2],x,y)

def remap_interval(val, input_interval_start, input_interval_end, output_interval_start, output_interval_end):
    """ Maps the input value that is in the interval [input_interval_start, input_interval_end]
        to the output interval [output_interval_start, output_interval_end].  The mapping
        is an affine one (i.e. output = input*c + b).

        Using math, we can find the formula necessary for remappying a value to another interval.
        Supposing we have a range [min,max] and we want to scale it to [a,b]

               (b-a)(x-min)
        f(x) = ------------  + a
                  max-min
    """
    output_interval = float(output_interval_end - output_interval_start)
    input_interval = float(input_interval_end-input_interval_start) #doesn't really matter which we float
    scaled_val = (output_interval*(val - input_interval_start)/(input_interval)) + output_interval_start
    return scaled_val

def create_image():
    """This function takes no inputs since we are generating a random image
    based on previously defined functions. It will generate separate functions
    and pixel values for each color channel and combine them in the final image.
    """
    red_func =  build_random_function(2,6)
    blue_func = build_random_function(2,4)
    green_func = build_random_function(2,5)

    img = Image.new("RGB", (350,350)) #we designate size of image here
    pixels = img.load()

    for x in range(350): #alter this value depending on size of image
        #Changing x value to within -1 and 1 to work with functions
        x_map = remap_interval(x,1,349,-1,1)
        for y in range(350):
            y_map = remap_interval(y,1,349,-1,1)
            #each set of code creates a random color map for the color specified
            red = evaluate_random_function(red_func,x_map,y_map)
            red_map = remap_interval(red,-1,1,0,255)

            green = evaluate_random_function(green_func,x_map,y_map)
            green_map = remap_interval(green,-1,1,0,255)

            blue = evaluate_random_function(blue_func,x_map,y_map)
            blue_map = remap_interval(blue,-1,1,0,255)
            #replaces values of the blank pixels in image with the generated color
            pixels[x,y] = (int(red_map),int(green_map),int(blue_map))
    img.show()

if __name__ == "__main__":
    print 'Builds random functions', build_random_function(2,5)
    print ' '
    print 'Evaluates random function', evaluate_random_function(['sin_pi', ['cos_pi', ['cos_pi', ['sin_pi', ['y']]]]],0.1,0.9)
    create_image()