# GGFW

Requires a working installation of ROOT (functionality to use a singularity container to run ROOT and avoid conflicts with ALICE software incoming).

To build the code and test that the library is working, run `sh build.sh`.

If you want to use the same folder organization software that I am using, run `sh setup.sh`. This will allow you to test the `AMPTExec.C` macro and the various histogram reader macros in the `draw` folder without changing too many file paths (though you will need to manually provide `AMPTExec.C` with the name and location of your AMPT input file when you add it to the chain on line 107. If you want to customize the name of the output file, it's line 106). 

Side note: at some point when this is more official and I have the time, I will change it so that there are instructions and more generic template draw files (same with `AMPTExec.C`). For now, I've just put the ones I am working on in their initial form so that you can get an idea of how to read out flow and pT correlations from the relevant containers - so take them with a grain of salt, they are not formatted super well!
