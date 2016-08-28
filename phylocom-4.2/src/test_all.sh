#/bin/tcsh

# test suite for phylocom

mkdir test >& /dev/null
cp ../phylo ../sample ../traits ../taxa .


echo '********* ./PHYLOCOM tests *************' > test/test.out
echo >> test/test.out

echo '========= phylo ===================' >> test/test.out
cat phylo >> test/test.out
echo >> test/test.out
echo 'cat phylo'

echo '========= sample ===================' >> test/test.out
cat sample >> test/test.out
echo >> test/test.out
echo 'cat sample'

echo '========= traits ===================' >> test/test.out
cat traits >> test/test.out
echo >> test/test.out
echo 'cat traits'

echo '========= ages ===================' >> test/test.out
cat ages >> test/test.out
echo >> test/test.out
echo 'cat ages'

echo '========= taxa ===================' >> test/test.out
cat taxa >> test/test.out
echo >> test/test.out
echo 'cat taxa'

echo '========= tree1 ===================' >> test/test.out
cat tree1 >> test/test.out
echo >> test/test.out
echo 'cat tree1'

echo '========= tree2 ===================' >> test/test.out
cat tree2 >> test/test.out
echo >> test/test.out
echo 'cat tree2'

echo '========= ./phylocom ===================' >> test/test.out
./phylocom >> test/test.out
echo >> test/test.out
echo 'phylocom'

echo '========= ./phylocom pd ===================' >> test/test.out
./phylocom pd >> test/test.out
echo >> test/test.out
echo 'phylocom pd'

echo '========= ./phylocom nodesig ===================' >> test/test.out
echo '  written to nodesig.nex' >> test/test.out
./phylocom nodesig > test/nodesig.out
echo >> test/test.out
echo 'phylocom nodesig'

echo '========= ./phylocom comstruct -m 0 ===================' >> test/test.out
./phylocom comstruct -m 0 >> test/test.out
echo >> test/test.out
echo 'phylocom comstruct -m 0'

echo '========= ./phylocom comstruct -m 1 ===================' >> test/test.out
./phylocom comstruct -m 1 >> test/test.out
echo >> test/test.out
echo 'phylocom comstruct -m 1'

echo '========= ./phylocom comstruct -m 2 ===================' >> test/test.out
./phylocom comstruct -m 2 >> test/test.out
echo >> test/test.out
echo 'phylocom comstruct -m 2'

echo '========= ./phylocom comstruct -m 3 -r 100 -w 50  ===================' >> test/test.out
./phylocom comstruct -m 3 -r 100 -w 50 >> test/test.out
echo >> test/test.out
echo 'phylocom comstruct -m 3'

echo '========= ./phylocom aot ===================' >> test/test.out
./phylocom aot >> test/test.out
echo >> test/test.out
echo 'phylocom aot'

echo '========= ./phylocom comdist ===================' >> test/test.out
./phylocom comdist >> test/test.out
echo >> test/test.out
echo 'phylocom comdist'

echo '========= ./phylocom phydist ===================' >> test/test.out
./phylocom phydist >> test/test.out
echo >> test/test.out
echo 'phylocom phydist'

echo '========= ./phylocom bladj ===================' >> test/test.out
./phylocom bladj >> test/test.out
echo >> test/test.out
echo 'phylocom bladj'

# echo '========= ./phylocom comnode ===================' >> test/test.out
# echo '  written to comnode.tre' >> test/test.out
# ./phylocom comnode > test/comnode.tre
# echo >> test/test.out
# echo 'phylocom comnode'

echo '========= ./phylocom sampleprune ===================' >> test/test.out
./phylocom sampleprune >> test/test.out
echo >> test/test.out
echo 'phylocom sampleprune'

echo '========= ./phylocom rndprune -p 10 -r 2===================' >> test/test.out
./phylocom rndprune -p 10 -r 2 >> test/test.out
echo >> test/test.out
echo 'phylocom rndprune'

echo '========= makenex ===================' >> test/test.out
./phylocom makenex > test/input.nex
echo '  written to input.nex'
echo >> test/test.out
echo 'phylocom makenex'

echo '========= ecovolve -h ===================' >> test/test.out
./ecovolve -h >> test/test.out
echo >> test/test.out
echo 'ecovolve -h'

echo '========= ecovolve ===================' >> test/test.out
./ecovolve >> test/test.out
echo '  also written to ecovolve1.new'
./ecovolve > test/ecovolve1.new
mv -f ecovolve.sample test/.
mv -f ecovolve.traits test/.
echo >> test/test.out
echo 'ecovolve'

echo '========= phylomatic -h ===================' >> test/test.out
./phylomatic -h >> test/test.out
echo >> test/test.out
echo 'phylomatic -h'

echo '========= phylomatic ===================' >> test/test.out
./phylomatic >> test/test.out
echo '  also written to phylomatic.new'
./phylomatic > test/phylomatic.new
echo >> test/test.out
echo 'phylomatic'

echo
echo 'DONE!'

rm -f sample traits phylo taxa

