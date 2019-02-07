library(data.table)
library(gtools)
library(ggplot2)
library(summarySE)

sorted_file_list <- function(path, pattern){
  f <- list.files(path=path, pattern=pattern, full.names=TRUE, recursive=TRUE)
  f2 <- gsub(".", " ", f, fixed=TRUE)
  f2_sort <- mixedsort(f2)
  f2_sort_rep <- gsub(" ", ".", f2_sort, fixed=TRUE)
  return(f2_sort_rep)
}

extract_info <- function(file, snp_pos, alt_base){
  f_spl <- gsub("[_/.]", " ", file)
  f_spl2 <- trimws(gsub(pattern='[A-Za-z]', "", f_spl, perl=TRUE))
  f_l <- as.numeric(strsplit(f_spl2, split=' ')[[1]])
  f_vec <- f_l[!sapply(f_l, is.na)][-(1:1)]
  
  f <- fread(file)
  f_rel <- f[snp_pos + 1, c(3,5,7,9)]
  
  obs <- f_rel[[alt_base]] / f_rel$cov
  rl <- f_vec[1]
  div <- f_vec[2]
  repl <- f_vec[3]
  sim <- f_vec[4]/100
  
  file_res <- c(rl, repl, div, sim, obs)
  return(file_res)
}

LM <- function(frame, div_val, rl_val, h){
  subframe <- subset(frame, div==div_val)
  model <- lm(subframe$sim ~ subframe$obs)
  slope <- model$coefficients[2][[1]]
  r <- summary(model)$adj.r.squared
  LM_text <- data.frame(lab=c(paste("r^2 == ", signif(r, digits=2), sep=''),
                              paste("slope == ", signif(slope, digits=2), sep='')),
                        sim=c(0.05, 0.05), obs=c(h, h-0.08), rl=rep(as.factor(rl_val),2),
                        div = rep(factor(paste(div_val,'% div.', sep='')), 2))
  return(LM_text)
}

#path_rl100 <- '~/path/to/sim/dirs/rl_100'
#path_rl1000 <- '~/path/to/sim/dirs/rl_1000'

dirs100 <- list.dirs(path=path_rl100, recursive=FALSE)
dirs1000 <- list.dirs(path=path_rl1000, recursive=FALSE)

rl100_list <- sorted_file_list(path=dirs100, pattern=".DOC5$")
rl1000_list <- sorted_file_list(path=dirs1000, pattern=".DOC5$")

# extract the info about every file
M100 <- t(sapply(rl100_list, extract_info, snp_pos=2000, alt_base='G'))
M1000 <- t(sapply(rl1000_list, extract_info, snp_pos=2000, alt_base='G'))

# turn to data.frame
AL_FREQ <- as.data.frame(rbind(M100, M1000))
AL_FREQ <- as.data.frame(M100)
names(AL_FREQ) <- c('rl', 'repl', 'div', 'sim', 'obs')

AL_FREQ <- subset(AL_FREQ, div!=0)
AL_FREQ$rl <- as.factor(AL_FREQ$rl)
AL_FREQ$div <- as.factor(AL_FREQ$div)

rl100_5 <- LM(subset(AL_FREQ, rl==100), div_val=5, rl_val=100, h=0.95)
rl100_10 <- LM(subset(AL_FREQ, rl==100), div_val=10, rl_val=100, h=0.95)
rl100_15 <- LM(subset(AL_FREQ, rl==100), div_val=15, rl_val=100, h=0.95)
rl100_20 <- LM(subset(AL_FREQ, rl==100), div_val=20, rl_val=100, h=0.95)
rl1000_5 <- LM(subset(AL_FREQ, rl==1000), div_val=5, rl_val=1000, h=0.81)
rl1000_10 <- LM(subset(AL_FREQ, rl==1000), div_val=10, rl_val=1000, h=0.81)
rl1000_15 <- LM(subset(AL_FREQ, rl==1000), div_val=15, rl_val=1000, h=0.81)
rl1000_20 <- LM(subset(AL_FREQ, rl==1000), div_val=20, rl_val=1000, h=0.81)

levels(AL_FREQ$div) <- list('5% div.' = '5', '10% div.' = '10', '15% div.' = '15', '20% div.' = '20')
AL_FREQ_SE <- summarySE(AL_FREQ, measurevar = 'obs', groupvars = c('rl', 'div', 'sim'))

plot.alfreq <- ggplot(AL_FREQ_SE, aes(x=sim, y=obs, color=rl)) +
  geom_smooth(method = "lm", size=0.5, alpha=0.2, se=FALSE) +
  ylim(c=0,1) +
  geom_point(data=AL_FREQ, alpha=0.3) +
  geom_point() +
  facet_grid(. ~ div) +
  geom_errorbar(aes(ymin=obs-se, ymax=obs+se), width=.05) +
  labs(color='Read length') + 
  ylab('observed') +
  xlab('simulated') +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  ggtitle('Allele frequency') +
  scale_color_manual(values=c('darkblue', 'darkred')) +
  geom_text(data=rl100_5, label=rl100_5$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl100_10, label=rl100_10$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl100_15, label=rl100_15$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl100_20, label=rl100_20$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl1000_5, label=rl1000_5$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl1000_10, label=rl1000_10$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl1000_15, label=rl1000_15$lab, hjust=0, parse=TRUE, size=3) +
  geom_text(data=rl1000_20, label=rl1000_20$lab, hjust=0, parse=TRUE, size=3) +
  theme_minimal() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=11))
plot.alfreq
ggsave(plot=plot.alfreq, width=20, height=9, unit="cm", filename='~/path/to/plot.pdf', device='pdf', dpi=300)









