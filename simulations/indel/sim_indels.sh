#! /bin/bash
# simulate reads for indel validation

chassis=$1
PGD='indel_valid.pgd'

build_pop_gen='../simulate_fork/build-population-genome.py'
read_sim='../simulate_fork/read_sim_pac.py'
te_seq='../te_seq/DOC5'


for j in 100 150 250 500 1000; do
	for i in {1..5}; do

		mkdir -p fastq
		
		python2.7 $build_pop_gen --chassis $chassis --pgd $PGD --output $PGD.fa


		GENOME=$(ls *.fa)
		e_rates=$(seq -w 0.00 0.02 0.32)

		for e in $e_rates; do
			epre=$(echo $e | cut -f2 -d'.')
			eprefix="0p$epre"'_'

			python2.7 $read_sim --pg $GENOME --read-length $j --haploid --fastq-prefix $eprefix --deletion-fraction 0.65 --error-rate $e
	
	        mv "$eprefix"'1.fastq' fastq
		done
		
		READS=$(ls fastq/*.fastq)
		for r in $READS; do
			echo $r
			deviaTE --input_fq $r --threads 4 --families 'DOC5' --library $te_seq
		done

		mkdir -p rl_${j}
		mkdir -p run_${i}
		mv fastq/*.DOC5 run_${i}
		mv run_${i} rl_${j}
		rm fastq/*

		
	done
done

echo "done"


