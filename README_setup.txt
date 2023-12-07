# On MAC OS M1 apple there is a little incompatibility issue where I need to give the location of stdio.h file
# I searched before using `locate stdio.h` line command
python setup.py build_ext -I/Library/Developer/CommandLineTools/SDKs/MacOSX12.3.sdk/usr/include build
