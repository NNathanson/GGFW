# GGFW

Requires a working installation of ROOT (functionality to use a singularity container to run ROOT and avoid conflicts with ALICE software incoming).

To build the code and test that the library is working, run `sh build.sh`.

If you want to use the same folder organization software that I am using, run `sh setup.sh`. This will allow you to test the `AMPTExec.C` macro and the various histogram reader macros in the `draw` folder without changing too many file paths (though you will need to provide it with the location of an AMPT input file). 
