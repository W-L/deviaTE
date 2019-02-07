#!/user/bin/env python
import os
import sys
import re
import time
import datetime
import random


class InsertRandomPopfreq:

    def __init__(self, genome_fasta="", te_fasta=""):
        self.genome_fasta = genome_fasta
        self.te_fasta = te_fasta

        
    def start(self, arg_dict):
        
        localtime = time.asctime( time.localtime(time.time()) )
        start_time = datetime.datetime.now().replace(microsecond=0)
        print ("\nCreating random TE insertion in Chasis genome started at\t"+str(localtime)+"\n")
        
        
        self.random_te_insertion(arg_dict)
        

        
        localtime = time.asctime( time.localtime(time.time()) )
        end_time = datetime.datetime.now().replace(microsecond=0)
        print ("\nCreating random TE insertion in Chasis genome completed at\t"+str(localtime))
        print ("Total time was taken \t"+str(end_time-start_time)+"\n")


    
    def random_te_insertion(self,arg_dict):
        
        
        
        insert_pos_list, seq_id, seq_length = self.get_chromosome_randon_insertion_pos(arg_dict)
        
        final_pop_freq_insertion_list = self.get_random_te_insertion(arg_dict)
      

        output_file = arg_dict["output"]
        
        fname, fext = os.path.splitext(output_file)
        temp_output_file = str(fname)+"_temp"
        
        ofh = open(temp_output_file,"w")
        
        #print (len(insert_pos_list))
        #print (len(final_pop_freq_insertion_list))
        
        if len(insert_pos_list)>1 and len(insert_pos_list)==len(final_pop_freq_insertion_list):
        
            index = 0
            
            for pos in insert_pos_list:
                final_pop_insert_str = final_pop_freq_insertion_list[index]
                to_print = "\t".join(map(str, [seq_id, pos, final_pop_insert_str]))
                ofh.write(str(to_print)+"\n")
                index = index+1
                
        ofh.close()
        
        if os.path.exists(temp_output_file) and os.path.isfile(temp_output_file):
            sort_cmd = "sort -k 1,1 -k 2,2n "+str(temp_output_file)+" > "+str(output_file)
            
            os.system(sort_cmd)
            os.remove(temp_output_file)
        
        
           
    
    
    
    def get_chromosome_randon_insertion_pos(self, arg_dict):
        flank_region = int(arg_dict["margin"])
        margin = int(arg_dict["flank_region"])
        
        te_fasta = arg_dict["te_fasta"]
        ref_fasta = arg_dict["ref_fasta"]
        
        te_seq_dict1, te_seq_dict2 = self.read_te_fasta(te_fasta)
        seq_id, seq_length, chr_seq = self.read_chesis_fasta(ref_fasta)
        
        te_count = len(te_seq_dict1)
        
                
        pop_size = int(arg_dict["pop_size"]) 
        pop_freq = arg_dict["pop_freq"] 
        insert_count = int(arg_dict["insert_count"]) 
        min_distance = int(arg_dict["min_distance"])
        
        
        min_dist_percent = 0
        min_dist_percent = int((min_distance*10)/100)
            
        pop_freq_list1 = pop_freq.split(",")
        pop_freq_list = []
        for freq in pop_freq_list1:
            freq = float(freq.lstrip().rstrip())
            if freq>0 and freq<=1:
                pop_freq_list.append(freq)
                
        freq_count = len(pop_freq_list)
        
        total_insertion = te_count * freq_count * insert_count
        
        new_min_dist = min_distance*total_insertion
        
        required_chr_seq_len = new_min_dist+flank_region+flank_region+(min_dist_percent*total_insertion)
        

        if (required_chr_seq_len>seq_length):
            sys.stderr.write("\tERROR: Chesis chromosome length is "+str(seq_length)+"bp and required length for insertion is "+str(required_chr_seq_len))
            sys.stderr.write("\n\n\tSuggestion:\n")
            sys.stderr.write("\tAdjust the --pop-freq, --insert-count and --min-distance parameters to meet the requirment\n\n")
            
            sys.exit(0)
            
            
        insert_pos_list = []
        
        i=flank_region
        while len(insert_pos_list)<=total_insertion:
            
            
            rand_min_distance = 0
            mu = min_distance # mean
            sigma = int(min_distance/min_dist_percent) # standard deviation
            
            
            start_pos = min_distance-sigma
            end_pos = min_distance+sigma
            #rn = random.lognormvariate(mu, sigma) # normal distribution random number generation
            
            #random.randrange(0, 101)
            rand_min_distance = random.randrange(start_pos, end_pos)
            
            #print (str(i), rand_min_distance, min_distance)
            #start_pos =
            #random.seed(i)
            #random.randint(a, b)
            
            insert_pos_list.append(i)
            i=i+rand_min_distance
            if len(insert_pos_list)>=total_insertion:
                break


        random.shuffle(insert_pos_list)
        random.shuffle(insert_pos_list)
        
        return insert_pos_list,seq_id, seq_length
    
    
    
    def get_random_te_insertion(self, arg_dict):
        
        te_fasta = arg_dict["te_fasta"]
        pop_size = int(arg_dict["pop_size"]) 
        pop_freq = arg_dict["pop_freq"] 
        insert_count = int(arg_dict["insert_count"]) 
        te_seq_dict1, te_seq_dict2 = self.read_te_fasta(te_fasta)
            
            
        pop_freq_list1 = pop_freq.split(",")
        pop_freq_list = []
        for freq in pop_freq_list1:
            freq = float(freq.lstrip().rstrip())
            if freq>0 and freq<=1:
                pop_freq_list.append(freq)
                    
        final_pop_freq_insertion_list = []
        for te in te_seq_dict1.keys():
            
            for freq in pop_freq_list:
                
                for ct in range(0,insert_count):
                    pop_to_insert = 0
                    pop_to_insert = int(float(freq)*int(pop_size))
                    
                    final_pop_insert_str = self.get_pop_rand_insertion(pop_size,te, pop_to_insert)

                    final_pop_freq_insertion_list.append(final_pop_insert_str)
        
        
        random.shuffle(final_pop_freq_insertion_list)
        random.shuffle(final_pop_freq_insertion_list)
        
        return final_pop_freq_insertion_list

    ## Step1 : read chesis chromosome
    def read_chesis_fasta(self,fasta_file):
        
        from fastaIO import FastaReader
        
        seq_dict = {}

        ct=1
        
        handle = FastaReader(fasta_file)
        
        seq_id=""
        seq_length=0
        chr_seq=""
        for seq_id,seq in handle:
            seq_length = len(seq)
            chr_seq=seq

        
        handle.close()
        
        return seq_id, seq_length, chr_seq
    
    
    ## Step2 : read TE fasta file
    def read_te_fasta(self,fasta_file):
        
        from fastaIO import FastaReader        
        seq_dict1 = {}
        seq_dict2 = {}

        ct=1
        handle = FastaReader(fasta_file)
        for seq_id,seq in handle:
            seq_dict1[ct] = seq_id
            seq_dict2[ct] = seq
            ct = ct+1
        
        handle.close()
        
        return seq_dict1, seq_dict2
    
    
    
    ## Step3: create random insertion population string number of column is equal to pop-size
    def get_pop_rand_insertion(self,pop_size,te, pop_to_insert):
        pop_num_list = []
        pop_col_list = []
        
        rand_pop = []
        for p in range(0,pop_size):
            pop_num_list.append(p)
            pop_col_list.append(0)
        
        random.shuffle(pop_num_list)
        random.shuffle(pop_num_list)
        rand_pop = pop_num_list[0:pop_to_insert]

        for rp in rand_pop:
            pop_col_list[rp]=te

        
        final_pop_insert_str = "\t".join(map(str, pop_col_list))
        #final_pop_insert_str1 = str(te)+"\t"+str(pop_to_insert)+"\t"+final_pop_insert_str
        return final_pop_insert_str
                    
    
        
    def utility(self, args, script_path):
        import multiprocessing
        

        ref_fasta = args.ref_fasta
        te_fasta = args.te_fasta
        pop_size = args.pop_size
        pop_freq = args.pop_freq
        output = args.output
        insert_count = args.insert_count
        min_distance = args.min_distance
        cpu_count = args.cpu_count
        
        
        if ref_fasta!=None and ref_fasta!="" and os.path.exists(ref_fasta) and os.path.isfile(ref_fasta):
            pass
        else:
            sys.stderr.write("\nERROR: --ref-fasta file path is not valid or does not exists\n\n")
            sys.exit(0)
        
        if te_fasta!=None and te_fasta!="" and os.path.exists(te_fasta) and os.path.isfile(te_fasta):
            pass
        else:
            sys.stderr.write("\nERROR: --te-fasta file path is not valid or does not exists\n\n")
            sys.exit(0)
            
            
        if cpu_count==None or cpu_count=="":
            cpu_count = int(multiprocessing.cpu_count())-2
            
        
        
        arg_dict = {}
        
        arg_dict["ref_fasta"] = ref_fasta
        arg_dict["te_fasta"] = te_fasta
        arg_dict["pop_size"] = pop_size
        arg_dict["pop_freq"] = pop_freq
        arg_dict["output"] = output
        arg_dict["insert_count"] = insert_count
        arg_dict["min_distance"] = min_distance
        arg_dict["script_path"] = script_path
        arg_dict["cpu_count"] = cpu_count


        fname, fext = os.path.splitext(arg_dict["ref_fasta"])
        combined_ref_genome = str(fname)+"+TE.fasta"
        if not os.path.exists(combined_ref_genome):
            cmd = "cat "+str(arg_dict["ref_fasta"])+" "+str(arg_dict["te_fasta"])+" > "+str(combined_ref_genome)
            os.system(cmd)
        
        flank_region=5 # flanking is used to add to the min_distance
        margin = 1000 # to leave 1000 base pair from the beginning and end of the chromosome for TE insertion
        
        arg_dict["margin"] = margin
        arg_dict["flank_region"] = flank_region
        
        return arg_dict