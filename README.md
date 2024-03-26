# EWF

An efficient simulator for exact Wright-Fisher diffusion and diffusion bridge paths, accounting for a wide class of selective regimes (genic, diploid and arbitrary polynomial selection), and the presence/absence of mutation. Please consult the _UserManual.pdf_ for all details with regards to installing dependencies, installing the program, and calling it from within python.

*Dependencies*

EWF requires the following:

- g++ compiler 
- boost library (https://boost.org)
- python and pip
- CMake
- pybind11

*Installation*

To install, run the following instructions in terminal at the root directory of EWF:
`mkdir build`   
`cd build`   
`cmake ..`   
`cmake --build .`   
`cd ..`   
`pip install .`

*Calling in python*

Once installed, you can call EWF from python by including   
`import EWF_pybind`    
at the start of your python script. Please see the scripts in the `examples` directory for more details and use cases.

If you come across any bugs, or have any queries/comments/suggestions please do get in touch using the email address Jaromir.Sant@gmail.com!