#!/usr/bin/env python3
import sys
import os.path
import subprocess
import itertools
from collections import OrderedDict
HOME = os.path.expanduser('~')

Params=OrderedDict()
exe = 'HartreeFock.exe'
hwlist=[25]
elist=[8]
e3list = [8]
file_extension_2n = ".me2j.gz"
file_extension_3n = ".me3j.gz"

emax_nn = 8
e2max_nn = 16
emax_3n =  8
e2max_3n = 8
e3max_3n = 8

Params['is_MBPTEnergy']=True
Params['type_3n_file']="no2b"
Params['dynamic_reference']=True
Params['alpha']=0.7
Params['iter_method']="mbroyden"
Params['iter_n_history']=9

element_table = [
    'NA',
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]
ZNList=[(2,2),(8,8),(20,20),(28,28),(50,50),(50,82),(62,82),(70,82),(70,102),(82,112),(82,126)]

def SetHamil(hw, emax_nn, e2max_nn, emax_3n, e2max_3n, e3max_3n):
    path2 = HOME + '/MtxElmnt/2BME'
    path3 = HOME + '/MtxElmnt/3BME'
    file_nn_int = 'TwBME-HO_NN-only_N4LO_EMN500_srg2.00_hw'
    file_nn_int += str(hw)+'_emax'+str(emax_nn)+'_e2max'+str(e2max_nn)+file_extension_2n
    #file_3n_int = "ThBME_srg2.00_N4LO_EMN500_ChEFT_N2LO_cD-1.80cE-0.31_LNL2_IS_hw"
    file_3n_int = "NO2B_ThBME_srg2.00_N4LO_EMN500_ChEFT_N2LO_cD-1.80cE-0.31_LNL2_IS_hw"
    file_3n_int += str(hw)+"_ms"+str(emax_3n)+"_"+str(e2max_3n)+"_"+str(e3max_3n)+file_extension_3n
    file_3n_int = "none"
    return path2, file_nn_int, path3, file_3n_int

def gen_script(params, file_nn, file_3n, batch, machine):
    fbase = "HF_"+params["Nucl"]+"_hw"+str(params["hw"])+\
            "_emax"+str(params["emax"])+"_e3max"+str(params["e3max"])
    fbase += "_" + file_nn.split(file_extension_2n)[0]
    fbase += "_" + file_3n.split(file_extension_3n)[0]
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    summary_file = "Summary_"+fbase+".dat"
    params["summary_file"] = summary_file
    fsh = "run_" + fbase + ".sh"
    prt = ""
    time = "0-04:00"
    if(params["emax"] <=16): time = "0-04:00"
    if(params["emax"] <=14): time = "0-02:00"
    if(params["emax"] <=12): time = "0-01:30"
    if(params["emax"] <=10): time = "0-01:00"
    if(params["emax"] <= 8): time = "0-00:45"
    if(params["emax"] <= 6): time = "0-00:30"
    if(params["emax"] <= 4): time = "0-00:20"
    if(machine=="oak"):
        prt += "#!/bin/bash \n"
        prt += "#PBS -q oak \n"
        prt += "#PBS -l mem=128gb,nodes=1:ppn=32,walltime=288:00:00 \n"
        prt += "cd $PBS_O_WORKDIR\n"
    if(machine=="cedar"):
        prt = "#!/bin/bash\n"
        prt += "#SBATCH --account="+account+"\n"
        prt += "#SBATCH --nodes=1\n"
        prt += "#SBATCH --ntasks=1\n"
        prt += "#SBATCH --cpus-per-task=32\n"
        prt += "#SBATCH --mem=125G\n"
        prt += "#SBATCH --time="+time+"\n\n"

    prt += 'echo "start ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value, str)):
            prt += str(key) + '= "' + str(value) + '" \n'
            continue
        if(isinstance(value, list)):
            prt += str(key) + "= "
            for x in value[:-1]:
                prt += str(x) + ", "
            prt += str(value[-1]) + "\n"
            continue
        prt += str(key) + '=' + str(value) + '\n'
    prt += "&end\n"
    prt += "EOF\n"
    if(batch):
        prt += exe + " " + file_input + " > " + file_log + " 2>&1\n"
        prt += "rm " + file_input + "\n"
    if(not batch):
        prt += exe + " " + file_input + "\n"
        prt += "rm " + file_input + "\n"
    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh

def get_item(params, Op="hamil"):
    if(not os.path.exists(params["summary_file"])): return None
    f = open(params["summary_file"],"r")
    lines = f.readlines()
    f.close()
    for line in lines:
        if(line[0] == "#"): continue
        try:
            data = line.split()
            if(data[0] != Op): continue
            return (float(data[1]), float(data[2]), float(data[3]), float(data[4]))
        except:
            return None

def main(machinename=None):
    if(machinename==None):
        batch = False
        machine = 'local'
    if(machinename!=None):
        batch=True
        if(machinename.lower() == "local"):
            machine = 'local'
        if(machinename.lower() == "oak"):
            machine = "oak"
        if(machinename.lower() =="cedar"):
            machine = "cedar"

    for ZN, e, e3, hw in itertools.product(ZNList,elist,e3list,hwlist):
        path2, file_nn, path3, file_3n = SetHamil(hw,emax_nn,e2max_nn,emax_3n,e2max_3n,e3max_3n)
        file_nn_int = path2 + '/' + file_nn
        file_3n_int = 'none'
        if(file_3n != 'none'):
            file_3n_int = path3 + '/' + file_3n

        Params['Nucl'] = element_table[ZN[0]] + str(ZN[0]+ZN[1])
        Params['hw']=hw
        Params['emax'] = e
        Params['e3max'] = e3
        Params['emax_nn'] = emax_nn
        Params['e2max_nn'] = e2max_nn
        Params['e2max_nn'] = e2max_nn
        Params['emax_3n'] = emax_3n
        Params['e2max_3n'] = e2max_3n
        Params['e3max_3n'] = e3max_3n
        Params['int_nn_file'] = file_nn_int
        Params['int_3n_file'] = file_3n_int
        Params['optrs']="Tkin,Hamil_NN_file_0+0,Hamil_NNN_file_0+0"
        Params['files_nn'] = "none,"+file_nn_int+",none"
        Params['files_3n'] = "none,none,"+file_3n_int
        fsh = gen_script(Params, file_nn, file_3n, batch, machine)
        if(machine == 'local'):
            cmd = "./" + fsh
        if(machine=="oak"):
            cmd = "qsub " + fsh
        if(machine=="cedar"):
            cmd = "sbatch " + fsh
        subprocess.call(cmd,shell=True)

    # create csv file
    prt = ""
    for ZN, e, e3, hw in itertools.product(ZNList,elist,e3list,hwlist):
        Params['Nucl'] = element_table[ZN[0]] + str(ZN[0]+ZN[1])
        Params['hw']=hw
        Params['emax'] = e
        Params['e3max'] = e3
        path2, file_nn, path3, file_3n = SetHamil(hw,emax_nn,e2max_nn,emax_3n,e2max_3n,e3max_3n)
        file_nn_int = path2 + '/' + file_nn
        file_3n_int = 'none'
        if(file_3n != 'none'):
            file_3n_int = path3 + '/' + file_3n
        fsh = gen_script(Params, file_nn, file_3n, batch, machine)
        en = get_item(Params)
        rp = get_item(Params,"Rp2")
        if(en == None): continue
        if(file_3n == "none"):
            prt += "{0:s}, {1:3d}, {2:3d}, {3:3d}, {4:3d}, {5:6.2f}, {6:12.6f}, {7:12.6f}, {8:12.6f}, {9:12.6f}, {10:12.6f}, {11:12.6f}, {12:12.6f}, {13:12.6f}\n".format(\
                    element_table[ZN[0]], ZN[0], ZN[1], ZN[0]+ZN[1], e, hw, en[0], en[1], en[2], en[3], rp[0], rp[1], rp[2], rp[3])
        if(file_3n != "none"):
            prt += "{0:s}, {1:3d}, {2:3d}, {3:3d}, {4:3d}, {5:3d}, {6:6.2f}, {7:12.6f}, {8:12.6f}, {9:12.6f}, {10:12.6f}, {11:12.6f}, {12:12.6f}, {13:12.6f}, {14:12.6f}\n".format(\
                    element_table[ZN[0]], ZN[0], ZN[1], ZN[0]+ZN[1], e, e3, hw, en[0], en[1], en[2], en[3], rp[0], rp[1], rp[2], rp[3])
    f = open("data.csv", "w")
    f.write(prt)
    f.close()

if(__name__ == '__main__'):
    if(len(sys.argv) == 1):
        main()
    if(len(sys.argv) > 1):
        main(sys.argv[1])


