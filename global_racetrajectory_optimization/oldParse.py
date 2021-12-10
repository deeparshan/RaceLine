import subprocess
import os
import configparser

my_path = os.getcwd()

if __name__=='__main__':
    lineNo = 0
    fileIn = open(my_path + '/params/racecar5.ini')
    optionS = []
    for line in fileIn.readlines():
        if (len(line.split(':'))==2):
            lineNo = lineNo + 1
            print(lineNo, line.rstrip())
            optionS.append(line.rstrip())
    
    option = input('Would you like to change anything? (y/n)')
    
    if (option == 'y'):
        try:
            choice = int(input('Enter the line you want to change'))
        except:
            print('Enter a number')
        
        lineNo = 0
        try:
            print(optionS[choice-1].split(':'))
        except:
            print("Set a valid range")    
                
            
    fileIn = open(my_path + '/params/racecar5.ini')
    fileIn = open(my_path + '/params/racecar5.ini')
    optionS = []
    for line in fileIn.readlines():
        if (len(line.split(':'))==2):
            print(line.rstrip())
            optionS.append(line.rstrip())
    try:
    	option = (input('Would you like to change anything? (y/n)'))
    except:
    	print('Wrong Choice.')
    		
    if option == "y":
    	option2 = ''
    	if (option == 'y'):
        try:
            choice = int(input('Enter the line you want to change'))
        except:
            print('Enter a number')
        
        lineNo = 0
        try:
            print(optionS[choice-1].split(':'))
        except:
            print("Set a valid range")    
                
    		
    elif option == "n":
    	print("No Change")
    else:
    	print('Invalid Option')