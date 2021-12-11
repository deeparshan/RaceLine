import subprocess
import os
from parseOptions import *
from main_compare import *

my_path = os.getcwd()

def optionB():
    completeEdit = False
    try:
        option1 = input('Would you like to change the two config files? (y/n)?: ')
        if (option1 == 'y' or option1 == 'n'):
            pass
    except:
    	print('Wrong Choice. Please enter y or n')
    if (option1 == 'y'):
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
        
        completeEdit = False
        while not (completeEdit):
            print("EDIT CAR 2")
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
    main_compare()
    return True

if __name__=='__main__':
	pass
