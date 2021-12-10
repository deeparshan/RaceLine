import subprocess
import os

my_path = os.getcwd()

if __name__=='__main__':
	subprocess.call("python3 ./parseOptions.py", shell=True)
