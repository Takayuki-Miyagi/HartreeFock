#!/usr/bin/env python3
import sys
import os.path
import subprocess

scriptname="run_hf_mbpt.py"

def gen_script( machine):
    fbase = "HF_serial_jobs"
    fsh = "run_" + fbase + ".sh"
    prt = ""
    if(machine=="oak"):
        prt += "#!/bin/bash \n"
        prt += "#PBS -q oak \n"
        prt += "#PBS -l mem=128gb,nodes=1:ppn=32,walltime=288:00:00 \n"
        prt += "cd $PBS_O_WORKDIR\n"
    if(machine=="cedar"):
        header = "#!/bin/bash\n"
        header += "#SBATCH --account="+account+"\n"
        header += "#SBATCH --nodes=1\n"
        header += "#SBATCH --ntasks=1\n"
        header += "#SBATCH --cpus-per-task=1\n"
        header += "#SBATCH --mem=125G\n"
        header += "#SBATCH --time=3-00:00\n\n"
    prt += "./"+scriptname + " local"
    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh


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

    fsh = gen_script(machine)
    if(machine == 'local'):
        cmd = "./" + fsh
    if(machine=="oak"):
        cmd = "qsub " + fsh
    if(machine=="cedar"):
        cmd = "srun " + fsh
    subprocess.call(cmd,shell=True)

if(__name__ == '__main__'):
    if(len(sys.argv) == 1):
        main()
    if(len(sys.argv) > 1):
        main(sys.argv[1])

