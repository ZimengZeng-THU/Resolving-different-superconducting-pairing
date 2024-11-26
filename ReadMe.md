Here I will introduce how to use the codes in this directory. With these codes, you can calculate the different superconductivity pairing channels of real materials base on epw calculation. 

1. You need to obtain the input files of the code in this directory. These input files are generated based on the modified QE code calculation. You can refer to the readme in the example directory. You also can download our file samples from the link[https://cloud.tsinghua.edu.cn/d/168a78a095d54d5d8865/] and direct start step next step.

2. Extract the input file(.lamda_kkq, .bxsf, .lambda_FS, wfc333.dat) to this directory, and then you need to reset the input of cal_delta.py based on the prefix and nkf1, nkf2, and nkf3 in epw.in. Running cal_delta.py will obtain the delta directory, which contains all the information about the gap functions.

3. You need to reset the input of irrep_delta.py based on the prefix and nkf1, nkf2, nkf3, efermi_energy in epw.in and the order if eigenvector you want to plot. Then you can run irrep_delta.py to analyze the representation of the delta functions and output the frmsf file, and then use fermisuffer to read this file for drawing. As a sample, you can run irrtp.delta.py and output leading odd parity channel result in default setting.