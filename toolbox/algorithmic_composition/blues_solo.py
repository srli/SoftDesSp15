""" Synthesizes a blues solo algorithmically """

"""NOTES FROM TESTING:
1. "scons" should be "scons setup.py" in order to generate the setup file
2. Running cvlc requires you to ctrl-C to stop running the program, is it supposed to do that?
3. Didn't need to install portaudio19-dev to run the code apparently... it wouldn't install on my system
4. It wasn't clear what a lick was in terms of Python datatypes. Do we want to introduce tuples?
I think having a few print statements to show what licks is would be helpful
Clarify that the first element of tuple changes the note, the second element changes the length of the note
5. Doing modular math while adding notes avoids indexing errors. Could let them figure that out though
Could clarify that they could do it in one line with modular math
"""

from Nsound import *
import numpy as np
from numpy.random import choice, rand, random_sample

def add_note(out, instr, key_num, duration, bpm, volume):
    """ Adds a note from the given instrument to the specified stream

        out: the stream to add the note to
        instr: the instrument that should play the note
        key_num: the piano key number (A 440Hzz is 49)
        duration: the duration of the note in beats
        bpm: the tempo of the music
        volume: the volume of the note
	"""
    freq = (2.0**(1/12.0))**(key_num-49)*440.0
    stream = instr.play(duration*(60.0/bpm),freq)
    stream *= volume
    out << stream

# this controls the sample rate for the sound file you will generate
sampling_rate = 44100.0
Wavefile.setDefaults(sampling_rate, 16)

bass = GuitarBass(sampling_rate)	# use a guitar bass as the instrument
solo = AudioStream(sampling_rate, 1)[(1,0.5), (1,0.5), (1, 0.5), (1, 0.5)]

""" these are the piano key numbers for a 3 octave blues scale in A
	See: http://en.wikipedia.org/wiki/Blues_scale """
blues_scale = [25, 28, 30, 31, 32, 35, 37, 40, 42, 43, 44, 47, 49, 52, 54, 55, 56, 59, 61]
beats_per_minute = 45				# Let's make a slow blues solo

curr_note = 0
add_note(solo, bass, blues_scale[curr_note], 1.0, beats_per_minute, 1.0)

licks = [ [(1,0.3), (2,0.7), (1, 0.4), (-2, 0.6)], [(0,0.5), (2,0.5), (-4, 0.5), (1, 0.5)], [(3,0.4), (-2,0.7), (-1, 0.6), (1, 0.3)] ]
#first element of tuple determines the note, the second element determins the length of the note
swing = True
for i in range(4): #we're running the lick 4 times, this range simply repeats the code 4 times
    lick = choice(licks)
    for note in lick:
        curr_note += note[0] % (len(blues_scale)-1) #we're moving up a single note
        if swing:
            length = 1.1
            swing = False
        elif not swing:
            length = 0.9
            swing = True
        note_length = note[1]*length
        add_note(solo, bass, blues_scale[curr_note], note_length, beats_per_minute, 1.0)

backing_track = AudioStream(sampling_rate, 1)
Wavefile.read('backing.wav', backing_track)

m = Mixer()

solo *= 0.4             # adjust relative volumes to taste
backing_track *= 2.0

m.add(2.25, 0, solo)    # delay the solo to match up with backing track
m.add(0, 0, backing_track)

m.getStream(500.0) >> "slow_blues.wav"