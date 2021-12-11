# RaceLine
KHADKAD, Deeparshan Khadka, ECE487: Project: Race Line Simulator
This document describes the run process of the program.
-------------

Here is the flow as of September 2021

Assumptions:
This is run on Ubuntu 18.02.

Requiremet:
There are a few functions that need to be installed to be able to run the Python code. There is a requirement.txt file in the module that can be used to get all the external functions. To do so,

Step1: go to the global_racetrajectroy_optimization directory (From ./Raceline/):
> cd global_racetrajectroy_optimization

Step2: Install all the python functions in requirements.txt, which can be done by the command line input
> pip3 install -r ./requirements.txt

NOTE: If there are any more installation errors, such as:
matplotlib requires tkinter -> can be solved by sudo apt install python3-tk
Python.h requires quadprog -> can be solved by sudo apt install python3-dev

------

With all installation requirements met, the program is now ready to be run.

STEP 1: Make sure that you are in global_racetrajectory_optimization
> cd ~/path/Raceline/global_racetrajectory_optimization

STEP 2: Run myApplication.py
> python3 ./myApplication.py

STEP 3: You will get prompted to a menu with 4 options
Make your pick, to determine the flow of the program, where your options are:
1. Run a custom car
2. Compare 2 Custom Cars
3. Simulate Last Run Setting 
4. Exit
------------
1. Run a custom car

If you decide to use menu option #1, you will be prompted to another menu that lists all the main configuration setup and asks you to choose what you'd like to change, and depending on if your chosen configuration has more dependancies within it, either another menu for those data gets prompted, or you will see the current value of the configuration and you'd be allowed to change it. This gets iterated until the user is satisfied with the new configuration. Then the simulation gets run.

------------
2. Compare 2 Custom Cars

Similarly to #1, here, you will be prompted to the car configuration menu twice for each car. You could opt to change only one or neither. Then, these two simulations get run.

------------
3. Simulate Last Run Setting

This will just run a single car without having the user go through the configuration files.

------------
4. Exit

Since this is a recursive menu, the user will get directed here every time the other 3 options finish their job. So, to exit, you will have to use the menu and the correct option.

-------------
Now, if you chose anything other than #4, you should see the simulation being run and the iterative function being called. Once the program calculates the optimal line, it starts outputting Matlab plots, which do not get saved automatically. The reason for it being that there are a lot of graphs that just overwhelm you with data, and saving plots using MatPlotLib gives a static output whereas just outputting the plot lets the user zoom in and out, and save portions that they think is neccessary.

However, the tajectory csv files should be saved in /global_racetrajectory_optimization/outputs

-------------
Other Data:
The global_racetrajectory_optimization has the README supplied by TUMFTM if you'd like to go through the working of thier program.
