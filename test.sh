# $1: either "clean" or sequence of processes the parallel
#     version should be executed with.
#
# $2: optional string of arguments for mpi executor defined
#     by $3
#
# $3: optional different mpi executor than mpirun
#
# $4: optional root directory

if [ ! "$4" == "" ]; then
  cd $4
fi

DIR=tests

l_seq="1 2 4 8 16 32 64 128 256 512"
seed_seq="1560 1561 1562 1563 1564"

if [ "$1" == "clean" ]; then
  rm -r $DIR
  echo Cleaned tests
  exit
fi

if [ "$3" == "" ]; then
  exc="mpirun"
else
  exc=$3
fi

if [ ! -d $DIR ]; then
  echo No tests found
  echo Generate correct outputs from serial version
  mkdir $DIR
  for l in $l_seq; do
    for seed in $seed_seq; do
      ./percolate_ser $seed -l $l --pgm_file_path "$DIR/ser.$l.$seed.pgm"
    done
  done
  echo Correct outputs generated
fi

echo Generate data from parallel
for i in $1; do
  for l in $l_seq; do
    for seed in $seed_seq; do
      echo Testing with i: $i, l: $l, seed: $seed
      $exc $2 -n $i ./percolate_par $seed -l $l --pgm_file_path "$DIR/par.$l.$seed.pgm"
      if [ -f "$DIR/par.$l.$seed.pgm" ]; then
        if [ "$(diff $DIR/par.$l.$seed.pgm $DIR/ser.$l.$seed.pgm)" == "" ]; then
          echo SUCCESS
          rm $DIR/par.$l.$seed.pgm
        else
          echo FAILURE
          #exit
        fi
      else
        echo NO FILE GENERATED
        exit
      fi
    done
  done
done
