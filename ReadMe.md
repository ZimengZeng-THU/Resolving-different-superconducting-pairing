# **Use of documentation of different pairing channels calculation of Pb**
Here I will introduce how to use the codes in this directory. With these codes, We can calculate the different superconductivity pairing channels of Pb base on epw calculation. 
## **Directory structure**
```text
Project Directory
├── cal_delta.py         # main function：calculate different pairing channels
├── irrep_delta.py       #  determined and output the gap function of  leading odd parity pairing channel of Pb
├── epw_calculation/     # the input file of QE to generate input data of cal_delta.py and `irrep_delta.py`
│   ├── scf.in           # the input file of QE
│   ...
│   └── Readme.md        # Use of documentation of how to generate input data of cal_delta.py and irrep_delta.py
├── result/              # output data and visualized result of cal_delta.py and irrep_delta.py
│   ├── delta/           # the output of `cal_delta.py` and the result(`irrep_delta_print.png` `Pb-Rep-gap.png`) of `irrep_delta.py`
│   │   ├── dkki-x.txt    # the eigenvector of the xth i-parity eigenvector of interaction matrix
│   │   ... 
│   ├── irrep_delta_print.png # irrep of C2x and inversion operation base on leading odd parity pairing channel
│   └── Pb-gap-rep.png   # Visualized the leading odd parity channel eigenvectors, corresponding to A1u in Fig.4 of paper
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
## **Input data and output of `cal_delta.py` and `irrep_delta.py`**
### **1.`cal_delta.py`**
- **Input data**：` .lambda_kkq` and `.lambda_FS`
- **Output data**: `delta/` directory
### **1.`irrep_delta.py`**
- **Input data**：` .bxsf`, `.lambda_FS` and `wfc333.dat`
- **Output data**: `pbt3.frmsf`
## **Calculation Porcess**
### **1.Generated of download input data**
First we should get input data of `cal_delta.py` and `irrep_delta.py`:
```bash
electron phonon coupling data and density of state data: .lambda_kkq; 
electron structure near fermi_surface: .lambda_FS; 
electron structure data: .bxsf;  
wavefunction data: wfc333.dat 
```
#### **generate data**
We can generate these input data by epw calculation,the process of this calculation can refer to `/epw_calculation/Readme.md`.
#### **download data**
We also can download these data in this link[https://cloud.tsinghua.edu.cn/d/168a78a095d54d5d8865/].

>**Note that:** this link also include directory `delta/`, it is output of `cal_delta.py`.
### **2.Calculate different pairing channel of Pb**
run the `cal_delta.py` code

```bash
 ./cal_delta.py
```

this code output the `delta/` directory, the `delta/dkki-x.txt` saved the xth i-parity eigenvector of interaction matrix.

>**Note that:** if we download input data from link[https://cloud.tsinghua.edu.cn/d/168a78a095d54d5d8865/], we don't need to run this step before starting next start.
### **3.Determined and output the leading odd parity channel**
We can run `irrep_delta.py` to analyze the representation of the delta functions and output the `.frmsf` file

```bash
./irrep_delta.py
```

the code should print character of irrep of C2x and inversion operation base on a pairing channel. 

![](./result/irrep_delta_print.png)

Because wannier function does not maintain perfect symmetry, the character here are not rigorous equal to 1. But in fact, the symmetry of wannier function is already good enough, so these character are almost 1 or -1.

This code also output the electron structure and leading parity gap function to `pbt3.frmsf` file. We can use FermiSurfer[https://mitsuaki1987.github.io/fermisurfer/] to open the `.frmsf` file for drawing.
### **4.visualized the eigenvector by fermisurfer**
We need to use FermiSurfer to open the `pbt3.frmsf`, set:
```text
Bar color: BMR
Background: 0  0  0
Line color: 1  1  1
```
and keep all other settings at their default values.

This result has shown in the A1u irrep of Fig(4) of our paper.
![](./result/Pb-gap-rep.png)

