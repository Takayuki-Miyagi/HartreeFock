#!/usr/local/bin/python3
import sys
import os.path
import subprocess
from get_filename import get_nnfile, get_nnnfile, get_summary_file,\
        get_hamil_file
from Conf import Orbit
HOME = os.path.expanduser('~')
params = {}
exe = './HartreeFock.exe'
hw = 35
renorm = 'srg'
cut = '2.00'
f2path = HOME + '/Desktop/TTFtest'
emax_2nf = 6
e2max_2nf = 12
pot = 'N3LO_EM500'
txtbin_2n = 'bin'
fom_2n = 'myg'
twbmefile = get_nnfile(f2path, 'TwBME', pot, renorm, cut, hw, \
                       emax_2nf, e2max_2nf, txtbin_2n, fom_2n)
scfile2 = 'None'

f3path = HOME + '/MtxElmnt/3BME'
emax_3nf = 14
e2max_3nf = 14
e3max_3nf = 14
e3cut = 14
genuine_3bf = True
# params for genuine n2lo 3BF
### R. Roth choice ##
cd = '-0.20'
ce = '0.098'
lambda_local = '400'
txtbin_3n = 'bin'
genuine_3bf = True
thbmefile = get_nnnfile(f3path, 'ThBME', renorm, cut, cd, ce, lambda_local, \
                        e3max_3nf, hw, genuine_3bf, txtbin_3n)
thbmefile = 'None'
scfile3 = 'None'

pmass = 8
nmass = 8
mass = 16
ENO = False
HFloop = False
thbme = False
if(thbmefile != 'None'): thbme = True
sv_hf_rslt = True
vac = 'ref'
MBPT = True
fmt_hf_snt = 'txt'
emax = 6
e2max = 12
conv = 1.e-8
params['hw'] = hw
params['emax_3nf'] = emax_3nf
params['e2max_3nf'] = e2max_3nf
params['e3max_3nf'] = e3max_3nf
params['e3cut'] = e3cut
params['pmass'] = pmass
params['nmass'] = nmass
params['mass'] = mass
params['sv_hf_rslt'] = sv_hf_rslt
params['vac'] = vac
params['fmt_hf_snt'] = fmt_hf_snt
params['HFloop'] = HFloop
params['emax'] = emax
params['e2max'] = e2max
params['conv'] = conv

name = 'hf.in'
f = open(name, 'w')
f.write('&input' + '\n')
for key, value in params.items():
    if(isinstance(value, str)):
        f.write(str(key) + '= "' + str(value) + '" \n')
    else:
        f.write(str(key) + '=' + str(value) + '\n')
f.write('&end' + '\n')
f.close()
ob = Orbit()
core = ob.get_core(nmass, pmass)
no = ob.get_no(nmass, pmass)
if(ENO):
    no = ob.get_no(nmass, pmass, 'ENO')
name1 = 'core.in'
f = open(name1, 'w')
f.write(str(len(core)) + '\n')
f.write('! n, l, j, itz \n')
for c in core:
    f.write(str(c[0]) + ' ' + str(c[1]) + ' ' + str(c[2]) + ' ' + str(c[3]) + '\n')
f.close()
name2 = 'no.in'
f = open(name2, 'w')
f.write(str(len(no)) + '\n')
f.write('! n, l, j, itz, fraction \n')
for o, fr in no.items():
    f.write(str(o[0]) + ' ' + str(o[1]) + ' ' + str(o[2]) + ' ' + str(o[3]) + ' ' + str(fr) + '\n')
f.close()
cmd = str(exe) + ' ' + name + ' ' + name1 + ' ' + name2 + ' ' + \
        twbmefile + ' ' + thbmefile + ' ' + scfile2 + ' ' + scfile3
print(cmd)
subprocess.call(cmd, shell=True)
f = get_summary_file(pmass, nmass, mass, HFloop, pot, renorm, cut, genuine_3bf, hw, emax, e2max, thbme)
cmd = 'mv Summary.out ' + f
subprocess.call(cmd, shell=True)
f = get_hamil_file(pmass, nmass, mass, HFloop, pot, renorm, cut, genuine_3bf, \
        hw, emax, e2max, fmt_hf_snt, thbme)
cmd = 'mv Hamil.snt.' + str(fmt_hf_snt) + ' '  + f
subprocess.call(cmd, shell=True)


