import sys
from optionA import *
from optionB import *
from optionC import *

menu_options = {
	1: 'Run a Custom Car',
	2: 'Compare 2 Custom Cars',
	3: 'Simulate Last Setting',
	4: 'Exit',
}

def print_menu():
	for key in menu_options.keys():
		print(key, '--', menu_options[key])
		
def option1():
    optionA()
	
def option2():
	optionB()
	
def option3():
	optionC()
    
def option4():
    sys.exit("Thank you for using my application")
    
	
if __name__=='__main__':
	while(True):
		print_menu()
		option = ''
		try:
			option = int(input('Your choice: '))
		except:
			print('Wrong Choice. Plese enter a number... ')
			
		if option == 1:
			option1()
		elif option == 2:
			option2()
		elif option == 3:
			option3()
		elif option == 4:
			option4()
		else:
			print('Invalid Option')
			

