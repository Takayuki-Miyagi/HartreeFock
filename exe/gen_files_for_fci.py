#!/usr/bin/env python3
import sys, time
import os.path
import subprocess
import itertools
from collections import OrderedDict
HOME = os.path.expanduser('~')

Params=OrderedDict()
hwlist=[16]
elist=[4]
e_tr_list=[4]
e3list = [16]
file_extension_2n = ".me2j.gz"
file_extension_3n = ".stream.bin"
ZNList=[(2,2)] # list of systems you are interested in

emax_nn  = 14   # this is for input NN filename
e2max_nn = 28  # this is for input NN filename
emax_3n  = 18  # this is for input 3N filename
e2max_3n = 36  # this is for input 3N filename
e3max_3n = 18  # this is for input 3N filename

operators = "R2,R4,Eccentricity2"
Params['optrs']=operators
Params['type_3n_file']="no2b"

# operator from file
#Params['files_nn'] = "file1,file2"
#Params['jpt_of_files'] = "jpz1,jpz2" example "0+0,1+1" (scalar and GT)

account="XXXX"         # Do not change
n_omp_threads=48       # threads per node
mem = "125G"           # memory per node
walltime = "0-01:00"   # request time
#scheduler = 'terminal'
scheduler = 'local'
#scheduler = 'slurm'
#scheduler = 'pbs'
if(scheduler == 'terminal' or scheduler == 'local'):
    header = ""
elif(scheduler == 'pbs'):
    header  = "#!/bin/bash \n"
    header += "#PBS -q oak \n"
    header += "#PBS -j oe \n"
    header += "#PBS -l mem=250gb \n"
    header += "#PBS -l walltime=72:00:00 \n"
    if(not mpi): header += f"#PBS -l nodes=1:ppn={n_omp_threads}\n"
    if(mpi): header += f"#PBS -l nodes={n_nodes}:ppn={n_omp_threads}\n"
    header += "cd $PBS_O_WORKDIR\n"
elif(scheduler == 'slurm'):
    header =  "#!/bin/bash\n"
    header += f"#SBATCH --account={account}\n"
    header += f"#SBATCH --nodes={n_nodes}\n"
    header += f"#SBATCH --ntasks-per-node=1\n"
    header += f"#SBATCH --cpus-per-task={n_omp_threads}\n"
    header += f"#SBATCH --mem={mem}\n"
    header += f"#SBATCH --time={walltime}\n\n"

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

def SetHamil(hw, emax_nn, e2max_nn, emax_3n, e2max_3n, e3max_3n):
    path2 = f"{HOME}/me2j" # directory of NN file
    path3 = f"{HOME}/me3j" # directory of 3N file
    file_nn_int = f"TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw{hw}_emax{emax_nn}_e2max{e2max_nn}.me2j.gz"
    file_3n_int = f"NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw{hw}_ms{emax_3n}_{e2max_3n}_{e3max_3n}.stream.bin"
    return path2, file_nn_int, path3, file_3n_int

def gen_script(params, batch):
    exe1 = 'HartreeFock.exe'
    exe2 = 'BasisTransform.exe'
    if(scheduler == 'slurm'):
        exe1 = "srun " + exe1
        exe2 = "srun " + exe2
    fbase = f"HF_{os.path.splitext(Params['OpFileName'])[0]}"
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    params["summary_file"] = f"summary_{fbase}.dat"
    fsh = "run_" + fbase + ".sh"

    prt = header
    if(scheduler == 'local'): prt += f'export OMP_NUM_THREADS={n_omp_threads}\n'
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
    if(scheduler == 'terminal'):
        prt += exe1 + " " + file_input + "\n"
        prt += exe2 + " " + file_input + "\n"
        prt += "rm " + file_input + "\n"
    else:
        prt += exe1 + " " + file_input + " > " + file_log + " 2>&1\n"
        prt += exe2 + " " + file_input + " >> " + file_log + " 2>&1\n"
        prt += "rm " + file_input + "\n"

    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh

def main():
    for ZN, e, e_tr, e3, hw in itertools.product(ZNList,elist,e_tr_list,e3list,hwlist):
        path2, file_nn, path3, file_3n = SetHamil(hw,emax_nn,e2max_nn,emax_3n,e2max_3n,e3max_3n)
        file_nn_int = path2 + '/' + file_nn
        file_3n_int = 'none'
        if(file_3n != 'none'):
            file_3n_int = path3 + '/' + file_3n

        Params['Nucl'] = element_table[ZN[0]] + str(ZN[0]+ZN[1])
        Params['Core'] = "n0"
        Params['hw']=hw
        Params['beta_cm'] = hw
        Params['emax'] = e
        Params['emax_mbpt']=e_tr
        Params['e3max'] = e3
        Params['emax_nn'] = emax_nn
        Params['e2max_nn'] = e2max_nn
        Params['e2max_nn'] = e2max_nn
        Params['emax_3n'] = emax_3n
        Params['e2max_3n'] = e2max_3n
        Params['e3max_3n'] = e3max_3n
        Params['int_nn_file'] = file_nn_int
        Params['int_3n_file'] = file_3n_int
        Params['OpFileName'] = f"{Params['Nucl']}_EM1.8_2.0_hw{Params['hw']}_eNAT{Params['emax']}_E{Params['e3max']}_e{Params['emax_mbpt']}.kshell.snt"
        Params['TransFileName'] = f"HOToNAT_{Params['Nucl']}_EM1.8_2.0_hw{Params['hw']}_e{Params['emax']}_E{Params['e3max']}.snt"
        Params['dynamic_reference']=False
        Params['alpha']=0.7
        Params['iter_method']="mbroyden"
        Params['iter_n_history']=9
        Params['is_NAT']=True
        Params['is_Op_out']=True
        Params['is_MBPT']=False

        fsh = gen_script(Params, False)
        if(scheduler == 'terminal'): cmd = "./" + fsh
        elif(scheduler == 'local'): cmd = "nohup ./" + fsh + " &"
        elif(scheduler == "pbs"): cmd = "qsub " + fsh
        elif(scheduler == 'slurm'): cmd = "sbatch " + fsh
        subprocess.call(cmd,shell=True)
        time.sleep(1)

if(__name__ == '__main__'):
    main()
