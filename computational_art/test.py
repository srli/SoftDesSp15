""" TODO: Put your header comment here """

import random
from math import cos, sin, pi
from random import randint
from PIL import Image
#from scitools.std import movie


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
    func = ['x','y','cos_pi','sin_pi'] #,'prod','square']
    if max_depth == 1:
        return hello[randint(0,1)]
    else:
        block = func[randint(2,3)]
        if block == 'prod': #accouts for when a block requires two inputs
            return [block, build_random_function(min_depth-1, max_depth-1), build_random_function(min_depth-1, max_depth-1)]
        elif not block == 'prod':
           return [block, build_random_function(min_depth-1, max_depth-1)]

# def build_random_function(min_depth, max_depth):
#     """ Builds a random function of depth at least min_depth and depth
#         at most max_depth (see assignment writeup for definition of depth
#         in this context)

#         min_depth: the minimum depth of the random function
#         max_depth: the maximum depth of the random function
#         returns: the randomly generated function represented as a nested list
#                  (see assignment writeup for details on the representation of
#                  these functions)
#     """
#     blocks = ["prod", "sin_pi", "cos_pi", "foo", "foobar", "avg"]
#     set_depth = random.randint(min_depth, max_depth)
#     def chooseBlock(depth):
#         curr_depth = depth + 1;
#         if curr_depth > set_depth:
#             #return "t"
#             return ["x",  "y", "t"][random.randint(0,2)]
#         choice = random.choice(blocks)
#         if choice=="prod" or "avg":
#             return [choice, chooseBlock(curr_depth), chooseBlock(curr_depth), chooseBlock(curr_depth)]
#         elif choice=="sin_pi" or "cos_pi":
#             return [choice, chooseBlock(curr_depth), chooseBlock(curr_depth)]
#         else:
#             return [choice, chooseBlock(curr_depth)]
#     res = chooseBlock(0)
#     return res

def evaluate_random_function(f, x, y, t=0):
    """ Evaluate the random function f with inputs x,y
        Representation of the function f is defined in the assignment writeup

        f: the function to evaluate
        x: the value of x to be used to evaluate the function
        y: the value of y to be used to evaluate the function
        returns: the function value
    """
    if f[0]=="x":
        return x
    elif f[0]=="y":
        return y
    elif f[0]=="t":
        return t
    elif f[0]=="prod":
        return evaluate_random_function(f[1], x, y, t) * evaluate_random_function(f[2], x, y, t) * evaluate_random_function(f[3], x, y, t)
    elif f[0]=="avg":
        return (1/3.0) * (evaluate_random_function(f[1], x, y, t) + evaluate_random_function(f[2], x, y, t) + evaluate_random_function(f[3], x, y, t))
    elif f[0]=="sin_pi":
        return sin(pi * evaluate_random_function(f[1], x, y, t) * evaluate_random_function(f[2], x, y, t))
    elif f[0]=="cos_pi":
        return cos(pi * evaluate_random_function(f[1], x, y, t) * evaluate_random_function(f[2], x, y, t))
    # elif f[0]=="foo":
    #     return evaluate_random_function(f[1], x, y, t)*5*t
    # elif f[0]=="foobar":
    #     return evaluate_random_function(f[1], x, y, t)**2


def remap_interval(val, input_interval_start, input_interval_end, output_interval_start, output_interval_end):
    """ Given an input value in the interval [input_interval_start,
        input_interval_end], return an output value scaled to fall within
        the output interval [output_interval_start, output_interval_end].

        val: the value to remap
        input_interval_start: the start of the interval that contains all
                              possible values for val
        input_interval_end: the end of the interval that contains all possible
                            values for val
        output_interval_start: the start of the interval that contains all
                               possible output values
        output_inteval_end: the end of the interval that contains all possible
                            output values
        returns: the value remapped from the input to the output interval

        >>> remap_interval(0.5, 0, 1, 0, 10)
        5.0
        >>> remap_interval(5, 4, 6, 0, 2)
        1.0
        >>> remap_interval(6, 4, 6, 0, 2)
        2.0
        >>> remap_interval(4, 4, 6, 0, 2)
        0.0
        >>> remap_interval(5, 4, 6, 1, 2)
        1.5
    """
    return (output_interval_end - output_interval_start) * 1.0 * (val - input_interval_start)/(input_interval_end - input_interval_start) + output_interval_start


def color_map(val):
    """ Maps input value between -1 and 1 to an integer 0-255, suitable for
        use as an RGB color code.

        val: value to remap, must be a float in the interval [-1, 1]
        returns: integer in the interval [0,255]

        >>> color_map(-1.0)
        0
        >>> color_map(1.0)
        255
        >>> color_map(0.0)
        127
        >>> color_map(0.5)
        191
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
    red_function = build_random_function(7,9)
    green_function = build_random_function(7,9)
    blue_function = build_random_function(7,9)

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

    im.save(filename)

# def generate_art(x_size=350, y_size=350):
#     """ Generate computational art and save as an image file.

#         filename: string filename for image (should be .png)
#         x_size, y_size: optional args to set image dimensions (default: 350)
#     """
#     # Functions for red, green, and blue channels - where the magic happens!
#     red_function = build_random_function(3,6)
#     green_function = build_random_function(3,5)
#     blue_function = build_random_function(3,6)
#     # Create image and loop over all pixels
#     num_frames = 120
#     filename = "art.png"
#     # filename = "./art"+str(k).zfill(3)+".png"
#     progress = remap_interval(k, 0, num_frames, -1, 1)
#     print(progress)
#     im = Image.new("RGB", (x_size, y_size))
#     pixels = im.load()
#     for i in range(x_size):
#         for j in range(y_size):
#             x = remap_interval(i, 0, x_size, -1, 1)
#             y = remap_interval(j, 0, y_size, -1, 1)
#             pixels[i, j] = (
#                     color_map(evaluate_random_function(red_function, x, y, progress)),
#                     color_map(evaluate_random_function(green_function, x, y, progress)),
#                     color_map(evaluate_random_function(blue_function, x, y, progress))
#                     )

#     im.save(filename)


if __name__ == '__main__':
    # import doctest
    # doctest.testmod()
    generate_art("art.png")
    # movie("./moveframes/frame*.png", fps=30, output_file="./movie.gif")
    # Create some computational art!
    # TODO: Un-comment the generate_art function call after you
    #       implement remap_interval and evaluate_random_function
    #generate_art("myart.png")

    # Test that PIL is installed correctly
    # TODO: Comment or remove this function call after testing PIL install
    #test_image("noise.png")