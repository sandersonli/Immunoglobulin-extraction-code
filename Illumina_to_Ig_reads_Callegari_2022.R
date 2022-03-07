#Illumina_to_Ig_reads_2021.R

# takes a list of sample identities in the format donor-well, and the filepath of a directory containing the results of an Illumina  paired end sequencing run 
# finds reads beloning to immunoglobulin chain genes, and writes them out as five text files (mu, gamma, alpha, kappa, lambda) in fasta format.
# The resulting fasta files can be assembled into contigs corresponding to the full chain genes using, for example, CAP3.

# gross structure of programme is:
#   for each antibody
#     for each chain
#         finds the indices of the constant hits
#         finds the indices of the variable hits for immunoglobulin gene that has more than representation_threshold (e.g., 5) constant hits
#         adds these indices to a vector of indices
#     next chain until all 5 done
#     for chain writes out a text file of the hits in fasta format
#   next antibody

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::valid()
BiocManager::install("SummarizedExperiment")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
library(ShortRead)
library(Biostrings)
library("SummarizedExperiment")

# the following three filepaths and the antibody list need to be specified by the user
# the script is currently configured to read 2 pairs of fastq.gz files corresponding to B_cell_transcriptome_abc and B_cell_transriptome_xyz in a directory called antibody_data
# the corresponding fast files will be called B_cell_transcriptome_abc_R1.fastq.gz, B_cell_transcriptome_abc_R2.fastq.gz, B_cell_transcriptome_xyz_R1.fastq.gz, B_cell_transcriptome_xyz_R2.fastq.gz
# the script uses a .csv file of constant region sequence fragments to start the search process, here provided at /antibody_data/multi_ig_probes.txt
# the user must provide a filepath to which the .fasta of identified reads will be saved, here a directory called immunoglobulin_sequences
# finally in line 203, a suffix can be specified for adding to each of the fasta files generated, here set to todays_date.txt. 

probe_length <- 17
multi_igs <-
  read.csv(file = "/antibody_data/multi_ig_probes.txt")
fastq_directory_filepath <-
  "/antibody_data/"
fasta_writeout_directory_filepath <-
  "/immunoglobulin_sequences/"
antibody_list <-
  c(
    "B_cell_transcriptome_abc",
    "B_cell_transcriptome_xyz"
  )

representation_threshold <- 5
Ig_genes <- c("mu", "gamma", "alpha", "kappa", "lambda")
########################################this section just prepares constant region probes for the 5 chains. The VC probes are already antisense.
mu_VC_probe <-"AAAGTGATGGAGTCGGGAAGGAAGTCCTGTGCGAGGCAGC" #40 base mu probe for trial

gamma_VC_probe <- "GAACCGGTGACGGTGTCGTGGAACTCAGGCGC|CTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGT"
alpha_VC_probe <- "GAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGGGGATGC|ACACTGAGTGGCTCCTGGGGGAAGAAGCCCTGGACCAGG"

kappa_VC_probe <- "CATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACT|CCTCTGTTGTGTGCCTGCTGAATAACTTCTATCCCAGAGA"

lambda_VC_probe <-"GTGACAGTGGCCTGGAAGGCAGATGGCAGCCCCGTCAAGG|CACTCTGTTCCCGCCCTCCTCTGAGGAGCTCCAAGCCAAC|CCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGG|ACCCGGGAGCCGTGACAGTGGCCTGGAAGGCAGAT"


mu_constant_frags <- DNAStringSet(multi_igs$mu_seqs)
mu_rev_frags <- reverseComplement(mu_constant_frags)
mu_constant_RC_probe <- paste(mu_rev_frags, collapse = "|")
mu_constant_probe <-
  paste(mu_constant_frags,
        mu_rev_frags,
        collapse = "|",
        sep = "|")
gamma_constant_frags <- DNAStringSet(multi_igs$gamma_seqs[multi_igs$gamma_seqs != "gamma"])
gamma_rev_frags <- reverseComplement(gamma_constant_frags)
gamma_constant_RC_probe <- paste(gamma_rev_frags, collapse = "|")
gamma_constant_probe <- paste(gamma_constant_frags, gamma_rev_frags, collapse = "|", sep = "|")
alpha_constant_frags <-
  DNAStringSet(multi_igs$alpha_seqs[multi_igs$alpha_seqs != "alpha"])
alpha_rev_frags <- reverseComplement(alpha_constant_frags)
alpha_constant_RC_probe <- paste(alpha_rev_frags, collapse = "|")
alpha_constant_probe <-
  paste(alpha_constant_frags,
        alpha_rev_frags,
        collapse = "|",
        sep = "|")
kappa_constant_frags <-
  DNAStringSet(multi_igs$kappa_seqs[multi_igs$kappa_seqs != "kappa"])
kappa_rev_frags <- reverseComplement(kappa_constant_frags)
kappa_constant_RC_probe <- paste(kappa_rev_frags, collapse = "|")
kappa_constant_probe <-
  paste(kappa_constant_frags,
        kappa_rev_frags,
        collapse = "|",
        sep = "|")
lambda_constant_frags <-
  DNAStringSet(multi_igs$lambda_seqs[multi_igs$lambda_seqs != "lambda"])
lambda_rev_frags <- reverseComplement(lambda_constant_frags)
lambda_constant_RC_probe <- paste(lambda_rev_frags, collapse = "|")
lambda_constant_probe <-
  paste(lambda_constant_frags,
        lambda_rev_frags,
        collapse = "|",
        sep = "|")
######################################## for each antibody in the provided antibody list extract the Ig hits from each of the 5 chains
for (antibody_sample in antibody_list) {
  antibody_R1_fastq_filepath <-
    paste(fastq_directory_filepath,
          antibody_sample,
          "_R1.fastq.gz",
          sep = "")
  antibody_R2_fastq_filepath <-
    paste(fastq_directory_filepath,
          antibody_sample,
          "_R2.fastq.gz",
          sep = "")
  R1_read_library <- readFastq(dirPath = antibody_R1_fastq_filepath)
  R2_read_library <- readFastq(dirPath = antibody_R2_fastq_filepath)
  
  ############################# for each of the 5 chains, find the constant, then the variable hits 
  for (Ig_gene in Ig_genes) {  ####FOR Ig genes
    if (Ig_gene == "mu") {
      ig_probe <- mu_constant_probe
      VC_probe <- mu_VC_probe
    }
    if (Ig_gene == "gamma") {
      ig_probe <- gamma_constant_probe
      VC_probe <- gamma_VC_probe
    }
    if (Ig_gene == "alpha") {
      ig_probe <- alpha_constant_probe
      VC_probe <- alpha_VC_probe
    }
    if (Ig_gene == "kappa") {
      ig_probe <- kappa_constant_probe
      VC_probe <- kappa_VC_probe
    }
    if (Ig_gene == "lambda") {
      ig_probe <- lambda_constant_probe
      VC_probe <- lambda_VC_probe
    }
    
    #################################################################### first the constant region hits
    R1_constant_hit_indices <-  grep(ig_probe, sread(R1_read_library))
    R2_constant_hit_indices <-  grep(ig_probe, sread(R2_read_library))
    total_hits <- length(R1_constant_hit_indices) + length(R2_constant_hit_indices)
    rolling_index_vector <- c(R1_constant_hit_indices, R2_constant_hit_indices)
    
    ############################################ IF C HITS only look for variable region hits if there are a > representation_threshold constant hits
    if (total_hits > representation_threshold) {
      ig_probe <- VC_probe
      # first get the VC hits from R2 library, and corresponding R1_V1 probes
      R2_hit_indices <- grep(ig_probe, sread(R2_read_library))
      R1_reads <- sread(R1_read_library)[unique(R2_hit_indices)]
      non_polyN_read_indices <- grep("A{10,}|C{10,}|G{10,}|T{10,}|N{2,}", invert = TRUE, R1_reads) ### take only reads without polyN, and more than 25 bases for probes
      R1_reads <- R1_reads[non_polyN_read_indices]
      R1_reads_full <- nchar(R1_reads) > probe_length+5
      print("before if")
      if (sum(R1_reads_full) > 0) {
        print("after R1 if")
        print(R1_reads)
        R1_reads_left_ends <- unique(substring(R1_reads[R1_reads_full], 5, probe_length+4))
        probes_from_R1 <- reverseComplement(DNAStringSet(R1_reads_left_ends))
      }
      
      
      # then get the VC hits from R1 library, and corresponding R2_V1 probes
      R1_hit_indices <- grep(ig_probe, sread(R1_read_library))
      R2_reads <- sread(R2_read_library)[unique(R1_hit_indices)]
      non_polyN_read_indices <- grep("A{10,}|C{10,}|G{10,}|T{10,}|N{2,}", invert = TRUE, R2_reads)
      R2_reads <- R2_reads[non_polyN_read_indices]
      R2_reads_full <- nchar(R2_reads) > probe_length+5
      if (sum(R2_reads_full) > 0) {
        print("after R2 if")
        R2_reads_left_ends <- unique(substring(R2_reads[R2_reads_full], 5, probe_length+4))
        print(length(R2_reads_left_ends))
        print(length(unique(R2_reads_left_ends)))
        probes_from_R2 <- reverseComplement(DNAStringSet(R2_reads_left_ends))
      }
      
      rolling_index_vector <- c(rolling_index_vector, R2_hit_indices, R1_hit_indices)
      
      ############ now, if probes were found, we want to just repeat the whole thing
      if (sum(c(R1_reads_full, R2_reads_full)) > 0){
        probes <- c(probes_from_R1, probes_from_R2)
        ig_probe <- paste(probes, collapse = "|")
        R2_hit_indices <- grep(ig_probe, sread(R2_read_library))

        # then get the VC hits from R1 library, and corresponding R2_V1 probes
        R1_hit_indices <- grep(ig_probe, sread(R1_read_library))

        rolling_index_vector <- c(rolling_index_vector, R2_hit_indices, R1_hit_indices)
        
        
      }########## end of if probes were found
      
    }########### end of IF C HITS
    
    rolling_index_vector <- unique(rolling_index_vector)
    all_hits_R1 <- R1_read_library[rolling_index_vector]
    all_hits_R2 <- R2_read_library[rolling_index_vector]
    
    fasta_writeout_filepath <-
      paste(
        fasta_writeout_directory_filepath,
        antibody_sample,
        "_",
        Ig_gene,
        "todays_date.txt",
        sep = ""
      )
    writeFasta(all_hits_R1, file = fasta_writeout_filepath, mode = "w")
    writeFasta(all_hits_R2, file = fasta_writeout_filepath, mode = "a")
    
    
  }########### end of FOR Ig gene
  
  
} 
######################################## next antibody



