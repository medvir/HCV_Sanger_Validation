library(tidyverse)
library(plotROC)

data = read.csv("/Volumes/huber.michael/Diagnostics/experiments/HCV_Sanger_Validation/data/HCV_Sanger_validation.csv") %>%
    #filter(FREQ <= 50) %>% ### only loking at minority mutations ?????
    select(pos, sample, Sanger, FREQ)


all_data=data.frame()
for (t in c(0, 1, 2, 5, 10, 15, 20, 25, 30, 40)) {
 
    data_i = data %>%
           
        mutate(TP = ifelse(FREQ >= t & Sanger == TRUE, 1, 0)) %>%   ### true positive (TP) eqv. with hit, recall
        mutate(FP = ifelse(Sanger == TRUE & FREQ == 0, 1, 0)) %>%  ### false positive (FP) eqv. with false alarm, Type I error
        mutate(FN = ifelse(Sanger == FALSE & FREQ >= t, 1, 0)) %>% ### false negative (FN) eqv. with miss, Type II error
        mutate(TN = ifelse(Sanger == FALSE & FREQ < t, 1, 0)) %>%  ### true negative (TN) eqv. with correct rejection
        mutate(checksum = TP + FP + FN + TN) %>%
        
        group_by(sample) %>% ### by sample
        
        mutate(TP = sum(TP, na.rm=TRUE)) %>%
        mutate(FP = sum(FP, na.rm=TRUE)) %>%
        mutate(FN = sum(FN, na.rm=TRUE)) %>%
        mutate(TN = sum(TN, na.rm=TRUE)) %>%
        
        top_n(1, pos*FREQ) %>%
        
        mutate(TPR = TP/(TP+FN)) %>%       ### sensitivity or true positive rate (TPR)
        mutate(SPC = TN/(TN+FP)) %>% ### specificity (SPC) or true negative rate
        mutate(FPR = 1-SPC) %>%      ### fall-out or false positive rate (FPR)
        
        select(sample, TP, FP, FN, TN, TPR, SPC) %>%
        mutate(threshold = t) %>%
        ungroup()
    
    all_data = rbind(all_data, data_i)
}

p = all_data %>%
    gather(key, value, 6:7) %>%
    ggplot(aes(x = threshold, y = value, color = key)) +
    geom_line() +
    facet_wrap( ~ sample) +
    geom_vline(xintercept = 15, color ="darkgrey") 
p