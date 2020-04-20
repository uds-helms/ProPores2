#!/bin/bash
echo
echo "Python 3 is required if you want to run the PROPORES user interface. It is NOT required when running PROPORES on the command line. Do you wish to install Python 3 and the required Python 3 packages (pyyaml and tkinter) on your system?"
echo
echo
read -p  "Type y for yes or n for no." -n 1
echo
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then 
	sudo apt-get install -y python3
	echo
	sudo apt-get install -y python3-tk
	echo
	sudo apt-get install -y python3-pip
	echo
	pip3 install --upgrade pip --user
	echo
	pip3 install pyyaml --user
	echo
fi
