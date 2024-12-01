# **Use of documentation of modified epw calculation**
Before calculate different pairing channels, we need to get the following data:
- **Electron-phonon coupling data**
- **Density of states (DOS) data**
- **Wavefunction data**
- **Electron structure data**

To achieve this, we modified the source code of EPW to enable the output of the corresponding data. The executable file compiled from the modified epw source code has been included in the `epw_calculation.zip`.
## **Directory structure**
```text
Project Directory
├── epw_calculation.zip  # Compressed QE executable file `pw.x`, `ph.x`, `epw-gkkgkk.x`, `epw-wfcreduce.x`
├── pb_s.UPF             # pseudopotential file of pb
├── scf.in               # the input file of self-consistent calculation
├── nscf.in              # the input file of non-self-consistent calculation
├── ph.in                # the input file of phonon calculation
├── epw.in               # input file of electron phonon coupling extraction by modified epw calculation
├── epw-w.in             # input file of wave function extraction by modified epw calculation
├── pp.py                # a code to extract all of input data about phonon calculation of epw calculation
├── weeph.py               # a code to extract wave function data in epw output directory
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
## **Calculation Porcess**
### **1.Unzip QE executable file**
Firstly, you need to unzip the epw_calculation zip file, which will give you pw.x, ph.x, epw_gkkgkk.x, and epw-wfcreduce.x. You need to run these executable files in the following step.
### **2.Phonon calculation**
The phonon code from QE requires a ground-state self-consistent run

`./pw.x scf.in <scf.out>`

Let us compute the dynamical matrix, the phonon frequencies and the variations of the self-consistent potential using the ph.x code from QE

`./ph.x ph.in <ph.out>`

Now we need to copy the input data of epw calculation (`.dyn`, `.dvscf`, and `.phsave` files) which have been produced by QE inside the `save/` folder.


`./pp.py`

### **3.Compute the Kohn-Sham wavefunctions on a coarse Brillouin zone grid**
In preparation of the EPW run we need to compute the Kohn-Sham wavefunctions and eigenvalues on a coarse Brillouin zone grid. For this we perform first a scf calculation and then an nscf calculation.

`./pw.x nscf.in <nscf.out>`

### **4.Compute electron phonon coupling matrix**
We use the modified EPW executable file `epw-gkkgkk.x` to compute and output the electron-phonon coupling matrix element data.

`./epw-gkkgkk.x epw.in <epw-gkkgkk.out>`

This step generates numerous output files. Our primary focus is on the `.lambda_kkq`(electron phonon coupling data and density of state data) , `lambda_FS` file(electron structure near fermi surface), `.bxsf` file(electron structure data) these is input data of different pairing channel calculation. 
### **5.Compute electron wave functions in wannier basis**
We use the modified EPW executable file `epw-wfcreduce.x` to compute and output the electron wave functions in wannier basis.

`./epw-wfcreduce.x epw.in <epw-wfc.out>`

This step generates numerous output files. Our primary focus is on the `ephmatx` files located in the `OUT/` directory. 

These files store the wavefunction projections, which are distributed across multiple files (`ephmat1`, `ephmat2`, etc.) based on different pools used during the calculation.

>**Note that:** that process will crash after output the wave function, don't care that

### **6.Save wave functions in `wfc333.dat`**
to deal the wave function in different pool

`./weeph.py`

this code output `wfc333.dat` to save wavefunction data.
### **7.Move some files to last directory**

`cp pb.lambda_kkq pb.lambda_FS pb.bxsf wfc333.dat ../`