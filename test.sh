make
./percolate -l 288 -d 0.411 -p 3

if [ "$(diff map.pgm true.pgm)" == "" ]; then
  echo SUCCESS
else
  echo FAILURE
fi
