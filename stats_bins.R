library(tidyverse)
library(readr)

file = "/Volumes/huber.michael/Diagnostics/experiments/HCV_Sanger_Validation/data/HCV_Sanger_validation.csv"

data = read.csv(file) %>% select(pos, sample, Sanger, FREQ)

bin_size = 5
alpha = 0.05
z = qnorm(1-alpha/2)

all_data=data.frame()
for (l1 in seq(0, 100-bin_size, bin_size)) {
        
        l2 = l1 + bin_size
        data_i = data %>%
                filter(FREQ > l1 & FREQ <= l2) %>%
                
                mutate(TP = ifelse(Sanger == TRUE, 1, 0)) %>% ### true positive (TP) eqv. with hit, recall
                mutate(FN = ifelse(Sanger == FALSE, 1, 0)) %>% ### false negative (FN) eqv. with miss, Type II error
        
                mutate(TP = sum(TP, na.rm=TRUE)) %>%
                mutate(FN = sum(FN, na.rm=TRUE)) %>%
                mutate(N = n()) %>%
                mutate(bin = paste0(">", l1, "-", l2)) %>%
                
                top_n(1, sample/pos) %>%
                
                mutate(TPR = round(TP/(TP+FN), digits = 4)) %>% ### sensitivity or true positive rate (TPR)
                mutate(p = (TP+0.5*z*z)/(N+z*z)) %>%  ### modified Wald Test
                mutate(lwr = round(p-z*sqrt(p*(1-p)/(N+z*z)), digits = 4)) %>% ### modified Wald Test
                mutate(upr = round(p+z*sqrt(p*(1-p)/(N+z*z)), digits = 4)) %>% ### modified Wald Test
                mutate(upr = ifelse(upr > 1, 1, upr)) %>%
                
                select(bin, N, TP, FN, lwr, TPR, upr) %>%
                mutate(bin = factor(bin, paste0(">", seq(0, 100-bin_size, bin_size),"-",seq(bin_size, 100, bin_size)))) 
        
        all_data = rbind(all_data, data_i)
}

p = all_data %>%
        ggplot(aes(x = bin, y = TPR, label = N)) +
        geom_point() +
        geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2) +
        geom_text(data=subset(all_data, TPR >= .95), aes(label = N), vjust = -1.5, hjust = NA) +
        geom_text(data=subset(all_data, TPR >= .2 & TPR < 0.95 ), aes(label = N), vjust = -4, hjust = NA) +
        geom_text(data=subset(all_data, TPR < .2), aes(label = N), vjust = -2, hjust = NA) +
        scale_y_continuous(breaks=seq(0,1,0.1))
p

write_excel_csv(all_data, "/Volumes/huber.michael/Diagnostics/experiments/HCV_Sanger_Validation/stats.xls")
