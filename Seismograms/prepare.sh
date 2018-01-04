nsta=30  #number of stations
ncom=3   #number of components


#put headers on top of recordings
#E-W component
  for i in $(seq 1 $nsta) ; do
    cat $(printf "head1.txt obs/obs_S%03i_C1.a" $i) > $(printf "st%02i.ve" $i)
  done
#N-S component
  for i in $(seq 1 $nsta) ; do
    cat $(printf "head1.txt obs/obs_S%03i_C2.a" $i) > $(printf "st%02i.vn" $i)
  done
#U-D component
  for i in $(seq 1 $nsta) ; do
    cat $(printf "head1.txt obs/obs_S%03i_C3.a" $i) > $(printf "st%02i.vz" $i)
  done

#SPECIAL CASES WITH DIFFERENT HEADERS!
 cat head3.txt obs/obs_S017_C1.a > st17.ve
 cat head3.txt obs/obs_S017_C2.a > st17.vn
 cat head3.txt obs/obs_S017_C3.a > st17.vz
 cat head2.txt obs/obs_S016_C1.a > st16.ve
 cat head2.txt obs/obs_S016_C2.a > st16.vn
 cat head2.txt obs/obs_S016_C3.a > st16.vz

