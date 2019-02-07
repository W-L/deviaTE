#! /bin/bash

# simulate reads for allele freq simulation

chassis=$1
read_len=$2

build_pop_gen='../simulate_fork/build-population-genome.py'
read_sim='../simulate_fork/read_sim_pac.py'
te_seq='../te_seq/DOC5'
divergences=$(seq 5 5 20)

mkdir -p pgds

for j in $divergences; do
	

	for i in {1..3}; do

		mkdir -p pgds
		python3 generate_pgd_allele_freq.py pgd_header_af $j
		
		mkdir -p fastq

		
		PGDs=$(ls pgds/*.pgd)
		for p in $PGDs; do
			echo $p
			python2.7 $build_pop_gen --chassis $chassis --pgd $p --output $p.fa
		done

		
		GENOMES=$(ls pgds/*.fa)
		for g in $GENOMES; do
			echo $g
			python2.7 $read_sim --pg $g --read-length $read_len --haploid --fastq-prefix $g
	        cat pgds/*.fastq >$g.fastq
	        mv $g.fastq fastq
	        rm pgds/*.fastq
		done

		
		READS=$(ls fastq/*.fastq)
		for r in $READS; do
			echo $r
			deviaTE --input_fq $r --threads 4 --families 'DOC5' --library $te_seq
		done
		
		mkdir -p run_$i
		mkdir -p div_$j
		mv fastq/*.DOC5 run_$i
		mv run_$i div_$j

		rm pgds/*
		rm fastq/*
	    
	done

done

echo "done"

