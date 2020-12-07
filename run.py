import os

for i in range(1, 9):
    os.system('sbatch --exclusive script/mpi_basic_weak/path' + str(i) + '.sub')
    os.system('sbatch --exclusive script/mpi_basic_strong/path' + str(i) + '.sub')
    os.system('sbatch --exclusive script/mpi_weak/path' + str(i) + '.sub')
    os.system('sbatch --exclusive script/mpi_strong/path' + str(i) + '.sub')
    os.system('sbatch --exclusive script/weak/path' + str(i) + '.sub')
    os.system('sbatch --exclusive script/strong/path' + str(i) + '.sub')
