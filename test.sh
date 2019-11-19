rm out/map.pgm

make

mpirun --hostfile hostfile -n $2 \
  $1 -p 3 --pgm_file_path "out/map.pgm"

echo Testing $1

if [ -f out/map.pgm ]; then
  if [ "$(diff out/map.pgm out/true.pgm)" == "" ]; then
    echo SUCCESS
  else
    echo FAILURE
  fi
else
  echo NO FILE GENERATED
fi
