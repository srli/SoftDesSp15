""" TODO: Put your header comment here """

from random import randint
from PIL import Image
from math import *

def build_random_function(min_depth, max_depth):
    """ Builds a random function of depth at least min_depth and depth
        at most max_depth (see assignment writeup for definition of depth
        in this context)

        min_depth: the minimum depth of the random function
        max_depth: the maximum depth of the random function
        returns: the randomly generated function represented as a nested list
                 (see assignment writeup for details on the representation of
                 these functions)
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
    """ Evaluate the random function f with inputs x,y
        Representation of the function f is defined in the assignment writeup
        f: the function to evaluate
        x: the value of x to be used to evaluate the function
        y: the value of y to be used to evaluate the function
        returns: the function value
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
    """ Given an input value in the interval [input_interval_start,

        input_interval_end], return an output value scaled to fall within
        the output interval [output_interval_start, output_interval_end].
    """
    output_interval = float(output_interval_end - output_interval_start)
    input_interval = float(input_interval_end-input_interval_start) #doesn't really matter which we float
    scaled_val = (output_interval*(val - input_interval_start)/(input_interval)) + output_interval_start
    return scaled_val


def color_map(val):
    """ Maps input value between -1 and 1 to an integer 0-255, suitable for
        use as an RGB color code.

        val: value to remap, must be a float in the interval [-1, 1]
        returns: integer in the interval [0,255]
    """
    # NOTE: This relies on remap_interval, which you must provide
    color_code = remap_interval(val, -1, 1, 0, 255)
    return int(color_code)


def test_image(filename, x_size=350, y_size=350):
    """ Generate test image with random pixels and save as an image file.

        filename: string filename for image (should be .png)
        x_size, y_size: optional args to set image dimensions (default: 350)
    """
    # Create image and loop over all pixels
    im = Image.new("RGB", (x_size, y_size))
    pixels = im.load()
    for i in range(x_size):
        for j in range(y_size):
            x = remap_interval(i, 0, x_size, -1, 1)
            y = remap_interval(j, 0, y_size, -1, 1)
            pixels[i, j] = (random.randint(0, 255),  # Red channel
                            random.randint(0, 255),  # Green channel
                            random.randint(0, 255))  # Blue channel

    im.save(filename)


def generate_art(filename, x_size=350, y_size=350):
    """ Generate computational art and save as an image file.

        filename: string filename for image (should be .png)
        x_size, y_size: optional args to set image dimensions (default: 350)
    """
    # Functions for red, green, and blue channels - where the magic happens!

    red_function =  build_random_function(2,20)
    blue_function = build_random_function(2,4)
    green_function = build_random_function(2,5)


    # red_function = ["x"]
    # green_function = ["y"]
    # blue_function = ["x"]

    # Create image and loop over all pixels
    im = Image.new("RGB", (x_size, y_size))
    pixels = im.load()
    for i in range(x_size):
        for j in range(y_size):
            x = remap_interval(i, 0, x_size, -1, 1)
            y = remap_interval(j, 0, y_size, -1, 1)
            pixels[i, j] = (
                    color_map(evaluate_random_function(red_function, x, y)),
                    color_map(evaluate_random_function(green_function, x, y)),
                    color_map(evaluate_random_function(blue_function, x, y))
                    )
    im.show()

generate_art("art.png")