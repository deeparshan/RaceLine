import subprocess
import os
from parseOptions import *
from main_globaltraj import *

my_path = os.getcwd()

def optionA():
    completeEdit = False
    while not (completeEdit):
        print("EDIT CAR 1")
        if (parseOptions()):
        	try:
        		option = input('Do you want to edit some more? (y/n): ')
        		if (option == 'y' or option == 'n'):
        			pass
        	except:
        		print('Wrong Choice. Please enter y or n')
        	if (option == 'n'):
        		completeEdit = True
    
    print("You are ready to run")
    main_globaltraj()
    return True


if __name__=='__main__':
	pass