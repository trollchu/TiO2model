make clean
make
mpirun -np 6  ./main inputs/TiO2_Rutile.dat output/data/densityandenergy.dat output/data/RDF.dat output/data/BAD.dat
make -C output/gnuplot/
