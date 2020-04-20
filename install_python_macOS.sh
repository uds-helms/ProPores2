#!/bin/bash
echo
echo "Python 3 is required if you want to run the PROPORES user interface. It is NOT required when running PROPORES on the command line. Do you wish to install Python 3 and the required Python 3 packages (pyyaml and tkinter) on your system? This will also install the package manager Homebrew if it is not installed yet."
echo
echo
read -p  "Type y for yes or n for no." -n 1
echo
echo

if [[ $REPLY =~ ^[Yy]$ ]]
then 
	which -s brew
	if [[ $? != 0 ]] ; then
		ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
		echo
	else
		brew update
		echo
	fi
	brew install python3
	echo
	pip3 install pyyaml --user
fi
