#!/usr/local/bin/python3
import sys
import os.path
import subprocess
HOME = os.path.expanduser('~')
class Orbit:
    def __init__(self):
        self.nlj_list = [{"n": 0, "l": 0, "j": 1},
                        {"n": 0, "l": 1, "j":  3},
                        {"n": 0, "l": 1, "j":  1},
                        {"n": 0, "l": 2, "j":  5},
                        {"n": 1, "l": 0, "j":  1},
                        {"n": 0, "l": 2, "j":  3},
                        {"n": 0, "l": 3, "j":  7},
                        {"n": 1, "l": 1, "j":  3},
                        {"n": 1, "l": 1, "j":  1},
                        {"n": 0, "l": 3, "j":  5},
                        {"n": 0, "l": 4, "j":  9},
                        {"n": 1, "l": 2, "j":  5},
                        {"n": 2, "l": 0, "j":  1},
                        {"n": 1, "l": 2, "j":  3},
                        {"n": 0, "l": 4, "j":  7},
                        {"n": 0, "l": 5, "j": 11},
                        {"n": 1, "l": 3, "j":  7},
                        {"n": 2, "l": 1, "j":  3},
                        {"n": 2, "l": 1, "j":  1},
                        {"n": 1, "l": 3, "j":  5},
                        {"n": 0, "l": 5, "j":  9},
                        {"n": 0, "l": 6, "j": 13},
                        {"n": 1, "l": 4, "j":  9},
                        {"n": 2, "l": 2, "j":  5},
                        {"n": 3, "l": 0, "j":  1},
                        {"n": 2, "l": 2, "j":  3},
                        {"n": 1, "l": 4, "j":  7},
                        {"n": 0, "l": 6, "j": 11},
                        ]

    def get_core(self,N,Z):
        core = []
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            if(num > N): continue
            nljtz = (n,l,j,1)
            core.append(nljtz)
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            if(num > Z): continue
            nljtz = (n,l,j,-1)
            core.append(nljtz)
        return core
    def get_valence(self,Nshell,Zshell):
        valence = []
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            if(2 * n + l != Nshell): continue
            nljtz = (n,l,j,1)
            valence.append(nljtz)
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            if(2 * n + l != Zshell): continue
            nljtz = (n,l,j,-1)
            valence.append(nljtz)
        return valence
    def get_no(self,N,Z,NO=None):
        no = {}
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            frac = 1
            if(num > N):
                if(NO != 'ENO'): break
                v = N - num + j + 1
                if(v > 0):
                    frac = float(v) / float(j+1)
                else:
                    frac = 0
                    break
            nljtz = (n,l,j,1)
            no[nljtz] = frac
        num = 0
        for nlj in self.nlj_list:
            n = nlj['n']
            l = nlj['l']
            j = nlj['j']
            num += j + 1
            frac = 1
            if(num > Z):
                if(NO != 'ENO'): break
                v = Z - num + j + 1
                if(v > 0):
                    frac = float(v) / float(j+1)
                else:
                    frac = 0
                    break
            nljtz = (n,l,j,-1)
            no[nljtz] = frac
        return no
params = {}
exe = './HartreeFock.exe'
hw = 35
renorm = 'srg'
cut = 2.0

f2path = './'
emax_2nf = 6
e2max_2nf = 12
pot = 'N3LO_EM500'
txtbin_2n = 'txt'
fom_2n = 'myg'
twbmefile = HOME + '/HF/TwBME-HO_NN-only_N3LO_EM500_srg2.00_hw35_emax6_e2max12.txt.myg'
twbmefile = 'None'
scfile2 = 'None'

f3path = './'
emax_3nf = 6
e2max_3nf = 6
e3max_3nf = 6
e3cut = 6
genuine_3bf = True
# params for genuine n2lo 3BF
### R. Roth choice ##
cd = -0.2
ce = 0.098
lambda_local = 400
txtbin_3n = 'txt'
thbmefile = HOME + '/HF/ThBME-srg2.00_cD-0.20cE0.098_lam400_e3max6_hw35_NNN-full.txt'
thbmefile = 'None'
scfile3 = 'None'

pmass = 8
nmass = 8
mass = 16
ENO = False
HFloop = True
HFbasis = True
NO2B = True
thbme = True
#thbme = False
sv_hf_rslt = True
vac = 'ref'
MBPT = True
fmt_hf_snt = 'bin'
emax = 3
e2max = 6
conv = 1.e-8
params['hw'] = hw
params['emax_2nf'] = emax_2nf
params['e2max_2nf'] = e2max_2nf
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
params['HFbasis'] = HFbasis
params['NO2B'] = NO2B
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
#valence = ob.get_valence(1, 1)
no = ob.get_no(nmass, pmass)
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
subprocess.call(cmd, shell=True)



