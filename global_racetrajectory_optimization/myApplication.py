import subprocess
import sys

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
    print('P1')
    subprocess.call("python3 ./optionA.py", shell=True)
	
def option2():
	subprocess.call("python3 ./parseOptions.py", shell=True)
	
def option3():
	print('OP3')
	subprocess.call("python3 ./optionC.py", shell=True)
	
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
			

