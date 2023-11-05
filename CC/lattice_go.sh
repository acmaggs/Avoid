#\rm -rf ./arun.*

for j in 1 2 3 4 5 6
do
for i in 1 2 3 4 5 6 
do
	mkdir brun.$i.$j
	cd brun.$i.$j
	../lattice &
	cd ..
done
done
wait
cat ?run.*/blist.dat > ./blist.dat
cat ?run.*/list.dat > ./list.dat
