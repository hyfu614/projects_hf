mpiexec -f mpd.hosts.265 -np 265 ./test_IRM_dsmh -F results/ -R test1_deg1 -H 50 -M 50 -G 264 -N 500 -D 1 1>tmp.out 2>tmp.err &
