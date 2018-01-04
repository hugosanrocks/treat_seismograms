nsta=17

  for i in $(seq 1 $nsta) ; do
#    cp $(printf "obs_S%03i_C1.a2 obs_S%03i_C1.a" $i $i)
    cp $(printf "obs_S%03i_C* ~/MEGA/INV3DKIN/run/jap/dat/" $i)
#    cp $(printf "obs_S%03i_C2.a2 obs_S%03i_C2.a" $i $i)
#    mv $(printf "obs_S%03i_C1.a ~/MEGA/INV3DKIN/run/jap/dat/" $i)
#    cp $(printf "obs_S%03i_C3.a2 obs_S%03i_C3.a" $i $i)
#    mv $(printf "obs_S%03i_C1.a ~/MEGA/INV3DKIN/run/jap/dat/" $i)
  done

