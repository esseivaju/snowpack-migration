#!/bin/bash

# use  nohup *.sh > fout  to redirect output to named file

# do valgrind with ctest !!!!
#TOOL="valgrind --tool=callgrind --simulate-cache=yes"
#TOOL="/software/bin/valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --log-fd=2 "
#TOOL="time"
#TOOL=""
#${TOOL} ../../bin/snowpack -c input/io_res1exp.ini -e 1996-06-17T00:00

../../bin/snowpack -c input/io_res1exp.ini -e 1996-06-17T00:00
