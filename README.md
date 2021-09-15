# BigBoyTelescope
Something I did in Tavatuy in 2019

# Basically, it's a program that calculates position of stars and the Moon
There are function in Main() that are commented out
Those only work if there's an Arduino that sends its GPS location.

# TXT FILES
Txt files are "The Database" of this project. It's a very bad way of storing data
But I couldn't do any better in the time that was given to us.
It wasn't our main focus anyway, we only needed to make the telescope move with the information that was sent to it.

# How it works
The Huh.ino file uploaded to an Arduino should be able to send GPS to .cpp program
In return it gets information about the angle of turning for two Stepper Motors (one X axis, one Y axis)

# How it works (2)
1.User points Telescope at the Moon

2.Program starts, calculates Moon's position

3.From GPS location it gives user a choice of stars, that are visible at the moment

4.User chooses the star

5.Program calculates star's coordinates relative to the Moon

6.Calculates the angle between the Moon and a chosen star

7.Sends information to Arduino

8.Arduino turns the Telescope using Stepper Motors

9.???

10.Profit, you now can see your desired star.
