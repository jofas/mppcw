rm tests/map.pgm

make

mpirun --hostfile hostfile -n $2 \
  $1 -p 3 --pgm_file_path "tests/map.pgm"

echo Testing $1

if [ -f tests/map.pgm ]; then
  if [ "$(diff tests/map.pgm tests/true.pgm)" == "" ]; then
    echo SUCCESS
  else
    echo FAILURE
  fi
else
  echo NO FILE GENERATED
fi
