# Shape Modeling and Geometry Processing projects

This repo contains the implementation of some algorithms in computational geometry.
## Projects Overview

### Project 1

Normals computation

### Project 2

Surface reconstruction

### Project 3

Differential Geonetry


## Installing Git and CMAKE
Before we can begin, you must have Git running, a distributed revision control system which you need to handin your assignments as well as keeping track of your code changes. We refer you to the online [Pro Git book](https://git-scm.com/book/en/v2) for more information. There you will also find [instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git]) on how to to install it. On windows we suggest using [git for windows](https://git-for-windows.github.io/).

CMake is the system libigl uses for cross-platform builds. If you are using Linux or macOS, I recommend installing it with a package manager instead of the CMake download page. E.g. on Debian/Ubuntu:
```
sudo apt-get install cmake
```
or with MacPorts on macOS:
```
sudo port install cmake.
```
On Windows you can download it from:
[https://cmake.org/download/](https://cmake.org/download/)



## Building Each Project
In the assignment repository you will find the different assignment directories 'assignmentX'. For now you only see the first one 'assignment1'. To compile the assignment code we will use the CMake building system.

Create a directory called build in the assignment directory, e.g.:
```
cd project1; mkdir build
```
Use CMake to generate the Makefiles needed for compilation inside the build/ directory:
```
cd build; cmake -DCMAKE_BUILD_TYPE=Release ../
```
On windows use the CMAKE gui with the buttons Configure and Generate.

Compile and run the executable, e.g.:
```
make && ./project1 <some mesh file>
```
Or use your favorite IDE.

Please compile your code in Release mode.

If you encounter any problems, please create an issue on the assignments repository on GitHub or ask the assistant in the exercise session.



