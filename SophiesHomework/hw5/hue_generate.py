# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 20:52:17 2014

@author: sophie
"""
from pattern.en import *
from pattern.web import *
import webcolors
from random import shuffle, randint
import Image
from quotes import quotes_pos, quotes_neg


def cleanup(stream):
    """Takes input from Twitter and Facebook. Removes common words so sentiment scores can be more accurate.
    Returns sentiment scores of each cleaned string of the input."""
    bad_words=['a','the','at','be','to','of','and','that','it','for','this','rt','is','in'] #common words we don't want
    clean_strings = []
    clean_sentiment = []
    i = 0
    while i < len(stream):
        string = stream[i]
        string_words = string.split(' ')
        clean_strings = ' '.join([j for j in string_words if j not in bad_words]) #only add word to clean_strings if it's not in bad words
        if sentiment(clean_strings)[0] != 0: #only add sentiment if the value isn't 0, people can't spell in twitter
            clean_sentiment.append(sentiment(clean_strings))
        i += 1
    return clean_sentiment

def twitter_call(hashtag):
    """Seaches twitter for the latest x tweets with the provided hashtag. Returns
    results in stream_tweet with no changes. Username search is possible, but oddly works on a unpredictable basis.
    This likely has something to do with twitter permissions"""
    t = Twitter()
    i = None
    stream_tweet = []
    for j in range(3):
        for tweet in t.search(hashtag, start=i, count=100): #searches twitter with input hashtag
            stream_tweet.append(tweet.text) #appends each tweet to stream_tweet
    if len(stream_tweet) < 1: #exits if twitter returns nothing due to twitter_call being unable to retrieve tweets sometimes
            raise Exception("This hashtag/twitter ID returns nothing. Please try another value.")
#    print stream_tweet #debugging
    twitter_sentiment = cleanup(stream_tweet) #runs sentiment analysis on cleaned tweets
    return twitter_sentiment
   
def facebook_call(facebook_license,name):
    """Searches facebook for lastest x posts and outputs to stream_facebook with no
    analysis. Ignores common posts by Facebook, i.e "Sophie Li has liked a photo" etc. """
    fb = Facebook(license=facebook_license, throttle=0.5)
    me = fb.profile(id=None) # (id, name, date, gender, locale, likes)-tuple
    stream_facebook = []
    for post in fb.search(me[0], type=NEWS, count=250):
        bad_status = [name, 'timeline', 'post', 'wall', 'status', 'photo'] #removes Facebook generated posts
        should_add = True
        for status in bad_status:
            if status in post.text: #only add the post to stream_facebook if it doesn't contain stuff in bad_status
                should_add = False
        if should_add:
            stream_facebook.append(post.text)
#    print stream_facebook #debugging
    facebook_sentiment = cleanup(stream_facebook) #runs sentiment analysis on cleaned facebook posts
    return facebook_sentiment
      
def sentiment_average(source_sentiment):
    """Takes list of all sentiments and averages all polarity and subjectivity scores"""
    positivity_total = 0
    subjectivity_total = 0
    i = 0
    while i < len(source_sentiment): #adds all polarity and subjectivity values together
        positivity_total = positivity_total + source_sentiment[i][0]
        subjectivity_total = subjectivity_total + source_sentiment[i][1]
        i += 1
    sentiment_average = [positivity_total/len(source_sentiment), subjectivity_total/len(source_sentiment)] #divides by total num values to find average
    return sentiment_average

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
    output_interval = output_interval_end - output_interval_start
    output_interval = float(output_interval)
    input_interval = input_interval_end-input_interval_start
    scaled_val = (output_interval*(val - input_interval_start)/(input_interval)) + output_interval_start
    return scaled_val

def closest_colour(requested_colour):
    """Generates RGB values that are closest to the ones in webcolor. Delivers answer by Euclidian distance in RGB space
    This function and get_colour_name were taken from fraxel on stackoverflow. 
    URL to answer is: http://stackoverflow.com/questions/9694165/convert-rgb-color-to-english-color-name-like-green"""
    min_colours = {}
    for key, name in webcolors.css3_hex_to_names.items():
        r_c, g_c, b_c = webcolors.hex_to_rgb(key)
        rd = (r_c - requested_colour[0]) ** 2
        gd = (g_c - requested_colour[1]) ** 2
        bd = (b_c - requested_colour[2]) ** 2
        min_colours[(rd + gd + bd)] = name
    return min_colours[min(min_colours.keys())]

def get_colour_name(requested_colour):
    """Outputs closest english color name to input RGB"""
    try:
        closest_name = actual_name = webcolors.rgb_to_name(requested_colour)
    except ValueError:
        closest_name = closest_colour(requested_colour)
        actual_name = None
    return actual_name, closest_name

def psycho_pass(facebook, name, twitter):
    """Pulls sentiment info in from twitter or Facebook. Takes all information and generates a color score"""
    if facebook == 'none': #allow user to input none if facebook or twitter doesn't apply
        facebook_sentiment = [(0,0.5)] #sets value to 0,0.5 because these are the middle values of sentiment analysis
    else:
        facebook_sentiment = facebook_call(facebook,name)
    
    if twitter == 'none':
        twitter_sentiment = [(0,0.5)]
    else:
        twitter_sentiment = twitter_call(twitter)
    
    total_sentiment = sentiment_average(facebook_sentiment + twitter_sentiment) #finds average sentiment across twitter and facebook   
    val1 = val2 = remap_interval(total_sentiment[0],-0.4,0.70,50,255) #maps 2 of RGB values to polarity
    val3 = remap_interval(total_sentiment[1],0.4,0.9,0,220) #since we're dealing with personal opinions, subjectivity is quite high
    #we have to modify scale to better characterize the smaller divisions in objectivity
    color_val = [int(val1), int(val2), int(val3)]
    shuffle(color_val) #shuffles values to allow generation of different hues
    color = (color_val[0], color_val[1], color_val[2])
    actual_name, closest_name = get_colour_name(color)
    #this 'quality' variable is a throwback to Psycho Pass, the inspiration for the project
    if total_sentiment[0] > 0.25:
        quality = 'clear'
    else:
        quality = 'steel'
    print ''
    print 'This individual\'s hue is:', quality, closest_name
    print ''
    print "Sentiment analysis shows that this individual or topic is", round(total_sentiment[0],2),"positive and", round(total_sentiment[1],2), "subjective"
    print ''
    img = Image.new("RGB", (350,350),color) #we designate size of image here
    img.show()
    return total_sentiment

def main():
    """This function runs when the user runs the script. Essentially pulls together all the functions for one cohesive package"""
    print "This script will look at your posts across social media and give you a value and color based on your general optimism"
    print "If you want to look at general Twitter sentiment about a topic, type in the desired hashtag"
    default = raw_input('Use preloaded settings? (Y/N): ').lower()
    if default == 'n':
        print 'Type none if entry does not apply. Entries are case sensitive. Do not put # in front of hashtags'
        name = raw_input('Enter your name as it appears on Facebook: ')
        facebook_license = raw_input('Please enter your Facebook license: ')
        twitter_username = raw_input('Please enter your Twitter username or a hashtag: ')
    else:
        print 'Default settings are Sophie Li\'s'
        name = 'Sophie Li'
        facebook_license = 'CAAEuAis8fUgBAN77AyE2NDJho0JJnP0ZBraJ5h7qmUwZBdtgJPeWb7tLSzg0NhxhDXBECuzmIKs11kujOnyfizyGxDgSlSB3OOiCMJYdCTTP7g9uarYGBCZC5t9HhDbHBdSo4HWGbOkIYX0MshCQT93BucAc0B3xKLUwZBmfsNoLuvLv3Pk6'
        twitter_username = 'me_pengineer'
    total_sentiment = psycho_pass(facebook_license, name, twitter_username)
    if total_sentiment[0] > 0.25: #prints quotes based on polarity. random thing I added for fun
        quote = randint(0,4)
        print 'In the words of', quotes_pos[quote][1],':'
        print quotes_pos[quote][0]
        print ''
        print 'Continue keeping bright things around.'
    else:
        quote = randint(0,4)
        print 'Remember the words of', quotes_neg[quote][1], 'who said:'
        print quotes_neg[quote][0]
        print ''
        print 'Look to happier things.'

if __name__ == '__main__':
    main()