library(tidyverse)
library(stringr)
library(seqinr)

path_to_files = "/Volumes/huber.michael/Diagnostics/experiments/HCV_Sanger_Validation/data/"

disambiguate <- function(nuc_list) {
        nuc_list = toupper(nuc_list[[1]])
        out = data.frame(pos = NULL, nt = NULL)
        for (i in 1:length(nuc_list)) {
                if (nuc_list[i] %in% c('A', 'C', 'G', 'T') ){
                        out = rbind(out, data.frame(pos = i, nt = nuc_list[i]))
                }
                else if (nuc_list[i] == "R" ){
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                }
                else if (nuc_list[i] == "Y" ){
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                }
                else if (nuc_list[i] == "M" ){
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                }
                else if (nuc_list[i] == "K" ){
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                }
                else if (nuc_list[i] == "S" ){
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                }
                else if (nuc_list[i] == "W" ){
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                }
                else if (nuc_list[i] == "H" ){
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                }
                else if (nuc_list[i] == "B" ){
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                }
                else if (nuc_list[i] == "V" ){
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                }
                else if (nuc_list[i] == "D" ){
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                }
                else if (nuc_list[i] == "N" ){
                        out = rbind(out, data.frame(pos = i, nt = "T"))
                        out = rbind(out, data.frame(pos = i, nt = "C"))
                        out = rbind(out, data.frame(pos = i, nt = "G"))
                        out = rbind(out, data.frame(pos = i, nt = "A"))
                }
        }
        return(out)
}


### Loop over all Sanger fasta files
files = list.files(path_to_files, pattern = ".fasta")

for (i in files) {
        name_i = gsub(".fasta", "", i)
        in_file = read.fasta(paste0(path_to_files, i))
        out_file = disambiguate(in_file) %>% mutate(Sanger = TRUE, sample = name_i)
        write.csv(out_file, paste0(path_to_files, name_i, ".csv"))
        out_fasta = out_file %>% group_by(pos) %>% top_n(1, nt)
        write.fasta(paste(out_fasta$nt, collapse = ""), name_i, paste0(path_to_files, name_i, "_disamb.fasta"))
        }

### Alignement of NGS reads to Sanger consensus





### Loop over all files
all_data=data.frame()
files = list.files(path_to_files, pattern = "lofreq.vcf")
for (i in files) {
        #if(i != "1000338554_lofreq.vcf") {next}
        name_i = gsub("_lofreq.vcf", "", i)
        
        vcf_file = paste0(path_to_files, name_i, "_lofreq.vcf")
        cons_file = paste0(path_to_files, name_i, "_NS5A_disamb.fasta")
        depth_file = paste0(path_to_files, name_i, ".depth")
        out_file = paste0(path_to_files, name_i, "_NS5A.csv")
        
        if (class(try(read.table(vcf_file))) == "try-error") {next} ### exit loop for this sample if vcf file is empty
        
        vcf_data = read.table(vcf_file, quote="\"") %>%
                rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, QUAL = V6, FILTER = V7, INFO = V8) %>%
                separate (INFO, c("DP", "AF", "SB", "DP4"), sep = ";", extra = "drop") %>%
                select(POS, REF, ALT, DP, AF) %>%
                mutate(DP = gsub("DP=" ,"", DP)) %>%
                mutate(AF = round(as.numeric(gsub("AF=" ,"", AF))*100,1))
        
        cons_data = data.frame(CONS = unlist(strsplit(readLines(cons_file)[-1], ""))) %>%
                mutate(POS = seq.int(nrow(.)))
        
        cov_data = read_delim(depth_file, "\t", col_names = FALSE, trim_ws = TRUE, col_types = "cii") %>%
                rename(POS = X2, COV = X3) %>%
                select(POS, COV)
        
        o = read.csv(out_file, header = TRUE) %>% select(pos, nt, Sanger)
        
        comb_data = full_join(cons_data, vcf_data, "POS") %>%
                full_join(cov_data, "POS") %>%
                mutate(REF = toupper(as.character(REF))) %>%
                mutate(ALT = toupper(as.character(ALT))) %>%
                mutate(CONS = toupper(as.character(CONS)))
        
        if (!(all(comb_data$CONS == comb_data$REF, na.rm = TRUE))) {next} ### exit loop for this sample if alignment not correct
        
        comb_data_2 = gather(comb_data, "ref_alt", "nt", 1:4, -REF, -POS, -AF, -COV, na.rm = TRUE) %>%
                select(ref_alt, POS, nt, AF, COV) %>%
                group_by(POS) %>%
                distinct(POS, nt, ref_alt, .keep_all = TRUE) %>%
                mutate(AF = ifelse(ref_alt == "CONS", 0, AF)) %>%
                mutate(sum_AF = sum(AF)) %>%
                mutate(FREQ = ifelse(ref_alt == "CONS", 100 - sum_AF, AF)) %>%
                rename(pos = POS) %>%
                select(-AF, -sum_AF, -ref_alt) %>%
                mutate(FREQ = ifelse(is.na(FREQ), 100, FREQ))
        
        master_table = full_join(o, comb_data_2, by = c("pos", "nt")) %>%
                mutate(Sanger = ifelse(is.na(Sanger), FALSE, Sanger)) %>%
                mutate(sample = name_i) %>%
                mutate(FREQ = ifelse(is.na(FREQ), 0, FREQ)) ### found by Sanger but not by NGS (false positive)
        
        
        all_data = rbind(all_data, master_table)
}

write.csv(all_data, paste0(path_to_files, "HCV_Sanger_validation.csv"))

### Plot
p1 = all_data %>%
        #filter(FREQ <= 50) %>%
        ggplot(aes(x=FREQ, y=Sanger, color = sample)) +
        geom_jitter() +
        theme(legend.position = "") +
        geom_vline(xintercept=15, color = "darkgrey")
p1

p2 = all_data %>%
        #filter(FREQ <= 50) %>%
        ggplot(aes(x=FREQ, y=Sanger, color = sample)) +
        geom_point(size = 0.5) +
        facet_wrap( ~ sample ) +
        theme(legend.position = "") +
        geom_vline(xintercept=15, color = "darkgrey")
p2

ggsave(paste0(path_to_files, "Figure_HCV_Sanger_validation_all.pdf"), p1, width = 30/2.54, height = 21/2.54)
ggsave(paste0(path_to_files, "Figure_HCV_Sanger_validation_facet.pdf"), p2, width = 30/2.54, height = 21/2.54)

