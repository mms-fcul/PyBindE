#!/usr/bin/python3
#python function to run bash commands
'''def bash_command(cmd):
    subprocess.Popen(cmd, shell=True, executable='/bin/bash')'''
    #function to use bash commands
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

#bash_command('a="Apples and oranges" && echo "${a/oranges/grapes}"')

#if the above doesn't work:
'''def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])


bash_command('g=/gromacs4/bin && echo $g')'''

#the way to do it with argparse:
'''parser = argparse.ArgumentParser()
parser.add_argument("echo", help="echo the string you use here")
args = parser.parse_args()
print(args.echo)
'''
#run_sas()
'''
gro="frame00000.gro"
tpr="sas.tpr"
ndx="pb.ndx"
sys="apolar_test"

run_sas(gro,tpr,ndx,sys)'''

'''
p = subprocess.Popen(['gmx', 'sasa','-f',gro,'-s',tpr,'-n',ndx,'-o',xvg],stdin=subprocess.PIPE)
p.communicate(b'Dimer\n')
p.wait()
'''

##def run_sas(gro, tpr, ndx, sysname):
    #sas_cmd = '{}/gmx sasa -f {} -s {} -n {} -o {}.xvg'.format(gmx20, gro, tpr, ndx, sysname)
    #subprocess.run(sas_cmd, shell=True)



'''def run_dynamics(mdp, gro, top, ndx, sysname):
    tpr_cmd = '{}/gmx grompp -f {} -c {} -p {} -n {} -o {}.tpr '\
              '-maxwarn 1000 -quiet'.format(md_configs.GroDIR, mdp, gro, top, ndx, sysname)
    sb.run(tpr_cmd, shell=True, stdout=md_configs.LOG, stderr=sb.STDOUT)
    mdrun_cmd = '{0}/gmx mdrun -nt {1} -s {2}.tpr -x {2}.xtc -c {2}.gro -e {2}.edr '\
                '-g {2}.log -o {2}.trr -rcon {3}'.format(md_configs.GroDIR, md_configs.nCPUs,
                                                         sysname, md_configs.rcon)
    sb.run(mdrun_cmd, shell=True, stdout=md_configs.LOG, stderr=sb.STDOUT)'''



