rm map.pgm

make
mpirun --hostfile hostfile --quiet -n 4 \
  percolate -l 4 -d 0.411 -p 3

if [ -f map.pgm ]; then
  if [ "$(diff map.pgm true.pgm)" == "" ]; then
    echo SUCCESS
  else
    echo FAILURE
  fi
else
  echo NO FILE GENERATED
fi
