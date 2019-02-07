#! /bin/bash

# simulate internal deletions for validation

chassis=$1

build_pop_gen='../simulate_fork/build-population-genome.py'
read_sim='../simulate_fork/read_sim_pac.py'
te_seq='../te_seq/DOC5'


for j in 80 100 175 250 500 1000; do
	mkdir -p pgds
	mkdir -p fastq

	python3 gen_pgd_int_del.py ./pgd_header

	# generate genomes from all PGDs
	PGDs=$(ls pgds/*.pgd)
	for p in $PGDs; do
		echo $p
		python2.7 $build_pop_gen --chassis $chassis --pgd $p --output $p.fa
	done
	
	# for each genome generate seuqencing reads, then merge the population sample
	GENOMES=$(ls pgds/*.fa)
	for g in $GENOMES; do
		echo $g
		python2.7 $read_sim --pg $g --read-length $j --haploid --fastq-prefix $g
	    cat pgds/*.fastq >$g.fastq
	    mv $g.fastq fastq
	    rm pgds/*.fastq
	done


	# run deviaTE on the sequencing reads
	READS=$(ls fastq/*.fastq)
	for r in $READS; do
		echo $r
		deviaTE --input_fq $r --threads 4 --families 'DOC5' --library $te_seq --no_freq_corr
	done

	mkdir rl_$j
	mv fastq/*.DOC5 rl_$j

	rm fastq/*
	rm pgds/*

done

echo "done"

