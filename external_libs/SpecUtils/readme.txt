SpecUtils is a meant for opening and manipulating spectrum files.
It was split off of InterSpec 20140820 in an attempt to be able to more easily
use this code in other projects.  But as a result of its evolution, it is
still a bit InterSpec centric, as well as rough and hacky; the CMake build 
script could also be improved.


To compile with support for Java use the following commands:
cd InterSpec/external_libs/SpecUtils
mkdir build
cd build
cmake -DSpecUtils_JAVA_SWIG=ON ..
make -j4

To then run the example Java executable, do:
cp ../swig/java_example/* .
javac -classpath .:jcommon-1.0.21.jar:jfreechart-1.0.17.jar:joda-time-2.9.jar *.java
java -Djava.library.path="./lib" -classpath .:jcommon-1.0.21.jar:jfreechart-1.0.17.jar:joda-time-2.9.jar Main

To compile with Python support, do:
cmake -DSpecUtils_PYTHON_BINDINGS=ON ..
