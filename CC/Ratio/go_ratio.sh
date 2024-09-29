(cd ..; make)
rm -f blist.dat
sleep 3
N=
iter=128
cp ../lattice .
for i in 1 2 3 4 5 6 7 8
do
    mkdir Run.$i
    cd Run.$i
    ../lattice $N $iter &
    cd ..
done
wait
cat */blist.dat >blist.dat
