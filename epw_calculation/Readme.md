Firstly, you need to unzip the epw_calculation zip file, which will give you pw.x, ph.x, epw_gkkgkk.x, and epw-wfcreduce.x. You need to run these executable files in the following order

to calculate phonon data
./pw.x scf.in <scf.out>
./ph.x ph.in <ph.out>

to get epw calculation input
./pp.py
./pw.x nscf.in <nscf.out>

to get electron phonon matrix file
./epw-gkkgkk.x epw.in <epw-gkkgkk.out>

to get wave function file
./epw-wfcreduce.x epw.in <epw-wfc.out>
that process will crash after output the wave function, don't care that

to deal the wave function in different pool
./weeph.py