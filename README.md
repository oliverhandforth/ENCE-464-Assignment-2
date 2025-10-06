ENCE461 Assignment 2 Template
=============================

See assignment instructions [here](doc/instructions/instructions.pdf)

Contents
--------
 - `doc/` - assignment instructions, lab notes, report template.
 - `reference/` - correct output for test comparison.
 - `main.c` - contains the command parsing logic and calls the poisson function.
 - `poisson.c` - basic template to work from. Write your solution here.
 - `test.sh` - automatic testing script.

Building
--------

Run `make` to build the poisson executable. Feel free to enhance the
makefile as needed, but the `poisson` executable must always be built through 
`make` .

Running
--------

Run the executable like so: `./poisson`

The following arguments are accepted:
* `-n`: the edge length of the cube.
* `-i`: the number of iterations.
* `-s`: the path/name of the file containing the source distribution
coordinates.
* `-t`: the number of threads to run. Once you develop your program, the
default
number of threads (in the absence of this argument) should be the most
optimal.

It also accepts the debug flag `-d` for verbose logging.

E.g. `./poisson -n 15 -i 100 -t 5 -s file.txt -d`

An arbitrary source distribution can be defined by listing the source points 
coordinates and their corresponding amplitudes in a text file in the following 
format:
```
x1,y1,z1,a1
x2,y2,z2,a2
....
```
Any undefined coordinates will have an amplitude of 0. 

You can find an 
example in `reference/asymmetrical-source-coordinates.txt`, which is used
in the test script.

In the absence of this argument, the default source distribution will be 
used, which is a single point with an amplitude of 1.0 at the center of 
the cube.


Testing
-------

Run `make test`.

It will automatically run your solution for different cube sizes and compare the
output against some correct reference files. **Do not edit these reference
files!**
