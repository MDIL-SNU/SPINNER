#!/bin/sh
#PBS -N test
#PBS -q batch
#PBS -l nodes=1:g13:ppn=28
#PBS -l walltime=200:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodefile
NPROC=`wc -l < $PBS_NODEFILE`

#############################################
################# change here ###############

inputfile='final.yaml'

csp_python='/home/kang1717/anaconda3/bin/python3.8'
training_python='python3.6'
amp2_python='python2.7'

num_of_core=$NPROC
mpi_command='mpirun'

#############################################
#############################################




##### READ INPUT #####

$csp_python src_final/read_input.py $inputfile
$csp_python src_final/print_num_gen.py $inputfile
ENCUT=`grep 'ENCUT' $inputfile\_read | tr -s ' ' | cut -d " " -f 3`
training_directory=`grep 'training_directory' $inputfile\_read | tr -s ' ' | cut -d " " -f 3`
input_dir=`grep 'input_dir' $inputfile\_read | tr -s ' ' | cut -d " " -f 2`
output_dir=`grep 'output_dir' $inputfile\_read | tr -s ' ' | cut -d " " -f 2`
total_iteration=`grep 'total_iteration' $inputfile\_read | tr -s ' ' | cut -d " " -f 3`
iteration_number=`grep 'number' $inputfile\_num_gen | tr -s ' ' | cut -d " " -f 2`
period=`grep 'period' $inputfile\_num_gen | tr -s ' ' | cut -d " " -f 2`
workdirectory=$PWD
##### MAIN ITERATION #####

for ((i=1; i<=iteration_number; i++))
do

# Crystal structure prediction
cd $workdirectory
rm -r src_final/__pycache__
cd src_final
$mpi_command -np $num_of_core $csp_python main.py $inputfile $i $period >& ../test.x
cd ..

$csp_python src_final/one_shot_monitoring.py $inputfile\_read $mpi_command $num_of_core $i

cd $output_dir
file=alert
if [ -a $file ]; then

   cd ../
  # DFT calculations
   # one shot
   $csp_python src_final/one_shot_calculation_final.py $inputfile\_read $mpi_command $num_of_core  $i
   # relax
   cd $workdirectory
   $csp_python src_final/sample_trainingset_final.py $inputfile\_read $mpi_command $i 

   cd $output_dir
   cd dft
   cat $pbs_nodefile > nodefile
   conf=./config.yaml
   src_final_path=`grep 'src_final_path' $conf | tr -s ' ' | cut -d " " -f 3`
   $amp2_python $src_final_path/main.py $conf nodefile $num_of_core $ENCUT >& stdout.x

   cp ERROR/*/POSCAR_* Submit
   rm -r ERROR/*
   conf=./config0.yaml
   src_final_path=`grep 'src_final_path' $conf | tr -s ' ' | cut -d " " -f 3`
   $amp2_python $src_final_path/main.py $conf nodefile $num_of_core $ENCUT >& stdout0.x

   # clean data and make str_list
   cd $workdirectory
   $csp_python src_final/clean_amp2_final.py $inputfile\_read $i 
   rm -r $output_dir/dft

   # Chage_directory name
   cd $workdirectory

   # NNP training
   cd $training_directory
   mkdir final$i
   cd    final$i

   cp $workdirectory/src_final/run.py .
   cp ../params* .
   cp $workdirectory/$output_dir/str_list$i .

   $csp_python $workdirectory/src_final/simple_nn_feature.py $workdirectory $output_dir\_iteration$i $i $inputfile
   cat str_list$i >> ../str_list_total
   cp ../str_list_total ./str_list
   $training_python run.py >& feature.out
   cp input.yaml input_feature.yaml

   $csp_python $workdirectory/src_final/simple_nn_train.py $workdirectory $output_dir\_iteration$i $i $inputfile
   cp ../scale_factor .
   cp ../pca .
   cp $workdirectory/$input_dir/potential ./potential_saved

   # end when converge criteria is met
   for ((j=1; j<=50; j++))
   do
   $training_python run.py >& training.out$j
   $csp_python $workdirectory/src_final/check_convergence_NNP.py $workdirectory/$inputfile

   file=nnp_done
   if [ -a $file ]; then
       break
   else
   cp potential_saved_iteration2001 potential_saved
   cp LOG LOG$j
   fi
   done

   cp potential_saved_iteration2001 $workdirectory/$input_dir/potential
   rm alert
fi
done
