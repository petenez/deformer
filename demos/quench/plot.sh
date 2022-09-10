for (( i = 0 ; i <= 10000 ; i += 1 )) ; do
  if [[ -e output-baboon-$i.rgb ]] && [[ ! -e baboon-quenched-$i.png ]] ; then
  	./pngwriter output-baboon-$i.rgb baboon-quenched-$i.png
  fi
done
