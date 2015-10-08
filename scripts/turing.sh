# @ job_name = smilei
# @ job_type = BLUEGENE
# Fichier sortie standard du travail
# @ output = $(job_name).$(jobid)
# Fichier erreur standard du travail
# @ error = $(output)
# Temps elapsed maximum demande
# @ wall_clock_limit = 10:00:00
# @ environment = JOBID=$(jobid)
# Taille bloc d'execution
# en nombre de noeuds de calcul (16 coeurs/noeud)
# @ bg_size = 2048
# @ queue


module load hdf5/mpi/1.8.14 python

mkdir OUT_$JOBID
cd OUT_$JOBID
runjob --ranks-per-node 16 --envs "OMP_NUM_THREADS=4" --envs "OMP_SCHEDULE=dynamic" --envs LD_LIBRARY_PATH=$BG_LD_LIBRARY_PATH PYTHONHOME=$BG_PYTHONHOME PYTHONPATH=$BG_PYTHONPATH --np 32768 : ../../src/smilei ../../tests/GC_sbs.in

