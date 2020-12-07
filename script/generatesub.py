numnode = 8

strongscriptdir = 'strong/'
weakscriptdir = 'weak/'
strongresultdir = 'results/strong/'
weakresultdir = 'results/weak/'

mpi_strong_script_dir = 'mpi_strong/'
mpi_weak_script_dir = 'mpi_weak/'
mpi_strongresultdir = 'results/mpi_strong/'
mpi_weakresultdir = 'results/mpi_weak/'

mydir = '$HOME/floyd-warshall'

# generate strong scaling .sub files

for i in range(0, numnode):
    filename = strongscriptdir + "path"+str(i+1)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ strongresultdir + 'path_'+str(i+1)+'n.out\n')
    f.write('#SBATCH -e '+ strongresultdir + 'path_'+str(i+1)+'n.err\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --ntasks=1\n')
    f.write('#SBATCH --tasks-per-node=1\n')
    f.write('#SBATCH --cpus-per-task=' +str(i+5)+'\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    f.write('\n')
    f.write('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n')
    f.write('cd ' + mydir + ' \n')
    f.write('./path.x -n 1800\n')
    f.write('\n')

# generate weak scaling .sub files for
for i in range(0, numnode):
    filename = weakscriptdir + "path"+str(i+1)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ weakresultdir + 'path_'+str(i+1)+'n.out\n')
    f.write('#SBATCH -e '+ weakresultdir + 'path_'+str(i+1)+'n.err\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --ntasks=1\n')
    f.write('#SBATCH --tasks-per-node=1 \n')
    f.write('#SBATCH --cpus-per-task=' +str(i+1)+'\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    f.write('\n')
    f.write('export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('module load openmpi-4.0.0 \n')
    f.write('cd ' + mydir + ' \n')
    f.write('./path.x -n ' + str((i+2)*300) + ' \n')
    f.write('\n')

#
# MPI-CANNON
#

# generate strong scaling .sub files for mpi
for i in range(1, numnode + 1):
    size = (i+1) * (i+1) + 1
    filename = mpi_strong_script_dir + "path"+ str(i)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ mpi_strongresultdir + 'path_'+str(i)+'n.out\n')
    f.write('#SBATCH -e ' + mpi_strongresultdir + 'path_' + str(i) + 'n.err\n')
    f.write('#SBATCH --ntasks=' + str(size) +'\n')
    f.write('#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    
    f.write('\n')

    f.write('cd ' + mydir + ' \n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('module load openmpi-4.0.0 \n')
    f.write('mpirun -n ' +  str(size)  + ' ./path-mpi-cannon.x -n 1800\n')
    f.write('\n')


# generate weak scaling .sub files for
for i in range(1, numnode + 1):
    size = (i+1) * (i+1) + 1

    filename = mpi_weak_script_dir + "path"+ str(i)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ mpi_weakresultdir + 'path_'+str(i)+'n.out\n')
    f.write('#SBATCH -e '+ mpi_weakresultdir + 'path_'+str(i)+'n.err\n')
    f.write('#SBATCH --ntasks=' + str(size) +'\n')
    f.write('#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    
    f.write('\n')

    f.write('cd ' + mydir + ' \n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('module load openmpi-4.0.0 \n')
    f.write('mpirun -n ' +  str(size) + ' ./path-mpi-cannon.x -n ' + str((i+1)*300) + '\n')
    f.write('\n')

#
# MPI-BASIC
#

mpi_basic_strong_script_dir = 'mpi_basic_strong/'
mpi_basic_weak_script_dir = 'mpi_basic_weak/'
mpi_basic_strongresultdir = 'results/mpi_basic_strong/'
mpi_basic_weakresultdir = 'results/mpi_basic_weak/'

# generate strong scaling .sub files for mpi
n = [5, 10, 18, 25, 36, 50, 72, 90]
tasks = [600, 900,1206,1500,1800,2100,2376, 2700]
for i in range(1, numnode+1):
    filename = mpi_basic_strong_script_dir + "path"+ str(i)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ mpi_basic_strongresultdir + 'path_'+str(i)+'n.out\n')
    f.write('#SBATCH -e '+ mpi_basic_strongresultdir + 'path_'+str(i)+'n.err\n')
    f.write('#SBATCH --ntasks=' + str(n[i-1])  +'\n')
    f.write('#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    
    f.write('\n')

    f.write('cd ' + mydir + ' \n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('module load openmpi-4.0.0 \n')
    f.write('mpirun -n ' +  str(n[i-1])  + ' ./path-mpi.x -n 1800\n')
    f.write('\n')


# generate weak scaling .sub files for
for i in range(1, numnode+1):
    filename = mpi_basic_weak_script_dir + "path"+ str(i)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J path\n')
    f.write('#SBATCH -o '+ mpi_basic_weakresultdir + 'path_'+str(i)+'n.out\n')
    f.write('#SBATCH -e '+ mpi_basic_weakresultdir + 'path_'+str(i)+'n.err\n')
    f.write('#SBATCH --ntasks=' + str(n[i-1])+'\n')
    f.write('#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    
    f.write('\n')

    f.write('cd ' + mydir + ' \n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('module load openmpi-4.0.0 \n')
    f.write('mpirun -n ' +  str(n[i-1]) + ' ./path-mpi.x -n ' + str(tasks[i-1]) + '\n')
    f.write('\n')
