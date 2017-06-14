
main:main_tio2.cpp ewald.cpp angular.cpp
	mpic++ main_tio2.cpp ewald.cpp angular.cpp -o main 

clean:
	rm -f main
	rm -f output/graph/*
	rm -f output/data/*

