make

mpirun --hostfile hostfile --quiet -n $2 \
  $1 -l 288 -d 0.411 -p 3 --pgm_file_path "out/map.pgm"

echo Testing $1

if [ -f out/map.pgm ]; then
  if [ "$(diff out/map.pgm out/true.pgm)" == "" ]; then
    echo SUCCESS
  else
    echo FAILURE
  fi
  rm out/map.pgm
else
  echo NO FILE GENERATED
fi
