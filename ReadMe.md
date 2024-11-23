Here I will introduce how to use the codes in this directory. With these codes, you can calculate the different pairing channels of real materials. 

1. you need to obtain the input files of the code in this directory. These input files are generated based on the modified QE code calculation. You can refer to the readme in the example directory or download our file samples from the link.

2. Extract the input file to this directory, and then you need to reset the input of cal_delta based on the prefix and nkf1, nkf2, and nkf3 in epw.in. Running cal_delta will obtain the delta directory, which contains all the information about the bandgap functions.

3. You need to run irrep_delta to analyze the representation of the delta functions and output the frmsf file, and then use fermisuffer to read this file for drawing