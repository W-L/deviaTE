#! /bin/bash

# builds the pop genome def. with python script first
# then runs the simulation 
# and analyses the results with deviaTE
# uses the chassis as argument

chassis=$1

build_pop_gen='../simulate_fork/build-population-genome.py'
read_sim='../simulate_fork/read_sim.py'
te_seq='../te_seq/DOC5'


for j in 100 150 250 500 1000; do
	echo "read length $j"
	
	mkdir -p pgds
	python3 generate_pgd_divergence.py

	for i in {1..5}; do
		
		
		PGDs=$(ls pgds/*.pgd)
		for p in $PGDs; do
			echo $p
			python2.7 $build_pop_gen --chassis $chassis --pgd $p --output $p.fa
			# break
		done
		echo ''
		

		GENOMES=$(ls pgds/*.fa)
		for g in $GENOMES; do
			echo $g
			python2.7 $read_sim --pg $g --read-length $j --haploid --fastq-prefix $g
			# break
		done
		echo ''

		
		READS=$(ls pgds/*.fastq)
		for r in $READS; do
			echo $r
			deviaTE --library $te_seq --threads 4 --families 'DOC5' --input_fq $r
			# break
		done
		echo ''

		mkdir -p rl_${j}
		mkdir -p run_${i}
		mv pgds/*.DOC5 run_${i}
		mv run_${i} rl_${j}
		
	done

done
echo "done"
