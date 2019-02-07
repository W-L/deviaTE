library(data.table)
library(gtools)
library(ggplot2)
library(summarySE)
# analyse divergence validation

file_list <- function(path, pattern){
  f <- list.files(path=path, pattern=pattern, full.names=TRUE, recursive=TRUE)
  f2 <- gsub(".", " ", f, fixed=TRUE)
  f2_sort <- mixedsort(f2)
  f2_sort_rep <- gsub(" ", ".", f2_sort, fixed=TRUE)
  return(f2_sort_rep)
}

fill_M <- function(files){
  f0 <- fread(files[1])
  f0_cov_len <- length(f0$cov)
  n_files <- length(files)
  M <- matrix(0, nrow = n_files, ncol = f0_cov_len)
  i <- 1
  
  for (file in files){
    f <- fread(file)
    M[i, ] <- f$cov
    i <- i +1
  }
  
  return(M)
}

MRD <- function(vec, pop_mean){
  # pop_mean = expected value
  dev_list <- lapply(vec, function(x){abs(x - pop_mean)/pop_mean})
  return(do.call(sum, dev_list)/length(dev_list))
}

snp_ratio <- function(files){
  perc_snp <- c()
  
  for (file in files){
    f <- fread(file)
    f_tab <- table(f$refsnp)
    total <- as.integer(f_tab[1])
    snps <- as.integer(f_tab[2])
    if (is.na(snps) == TRUE){
      snps <- 0
    }
    perc_snp <- c(perc_snp, (snps/total) * 100)
  }
  return(perc_snp)
}

pth <- '~/path/to/divergence/validation/'

rl100_list <- file_list(path=paste(pth, 'divergence/rl_100/', sep=''), pattern=".DOC5$")
rl150_list <- file_list(path=paste(pth, 'divergence/rl_150/', sep=''), pattern=".DOC5$")
rl250_list <- file_list(path=paste(pth, 'divergence/rl_250/', sep=''), pattern=".DOC5$")
rl500_list <- file_list(path=paste(pth, 'divergence/rl_500/', sep=''), pattern=".DOC5$")
rl1000_list <- file_list(path=paste(pth, 'divergence/rl_1000/', sep=''), pattern=".DOC5$")

rl100_list_indel <- file_list(path=paste(pth, 'divergence_indel/rl_100/', sep=''), pattern=".DOC5$")
rl150_list_indel <- file_list(path=paste(pth, 'divergence_indel/rl_150/', sep=''), pattern=".DOC5$")
rl250_list_indel <- file_list(path=paste(pth, 'divergence_indel/rl_250/', sep=''), pattern=".DOC5$")
rl500_list_indel <- file_list(path=paste(pth, 'divergence_indel/rl_500/', sep=''), pattern=".DOC5$")
rl1000_list_indel <- file_list(path=paste(pth, 'divergence_indel/rl_1000/', sep=''), pattern=".DOC5$")


# load the data into a matrix: row = coverage in one file, colum = pos
M100 <- fill_M(rl100_list)
M150 <- fill_M(rl150_list)
M250 <- fill_M(rl250_list)
M500 <- fill_M(rl500_list)
M1000 <- fill_M(rl1000_list)

M100_indel <- fill_M(rl100_list_indel)
M150_indel <- fill_M(rl150_list_indel)
M250_indel <- fill_M(rl250_list_indel)
M500_indel <- fill_M(rl500_list_indel)
M1000_indel <- fill_M(rl1000_list_indel)

# calc mean absolute deviation from the expected value
M100_MRDs <- apply(M100, 1, MRD, 100) * 100
M150_MRDs <- apply(M150, 1, MRD, 150) * 100
M250_MRDs <- apply(M250, 1, MRD, 250) * 100
M500_MRDs <- apply(M500, 1, MRD, 500) * 100
M1000_MRDs <- apply(M1000, 1, MRD, 1000) * 100

M100_MRDs_indel <- apply(M100_indel, 1, MRD, 100) * 100
M150_MRDs_indel <- apply(M150_indel, 1, MRD, 150) * 100
M250_MRDs_indel <- apply(M250_indel, 1, MRD, 250) * 100
M500_MRDs_indel <- apply(M500_indel, 1, MRD, 500) * 100
M1000_MRDs_indel <- apply(M1000_indel, 1, MRD, 1000) * 100

div_col <- rep(seq(0,30), 5)
indel_col <- rep(seq(0,32, by=2), 5)


MRD_div100 <- data.frame('div'=div_col, 'MRD'=M100_MRDs, 'read_len'=rep('100', length(div_col)), 'type'=rep('Mismatches', length(div_col)))
MRD_div150 <- data.frame('div'=div_col, 'MRD'=M150_MRDs, 'read_len'=rep('150', length(div_col)), 'type'=rep('Mismatches', length(div_col)))
MRD_div250 <- data.frame('div'=div_col, 'MRD'=M250_MRDs, 'read_len'=rep('250', length(div_col)), 'type'=rep('Mismatches', length(div_col)))
MRD_div500 <- data.frame('div'=div_col, 'MRD'=M500_MRDs, 'read_len'=rep('500', length(div_col)), 'type'=rep('Mismatches', length(div_col)))
MRD_div1000 <- data.frame('div'=div_col, 'MRD'=M1000_MRDs, 'read_len'=rep('1000', length(div_col)), 'type'=rep('Mismatches', length(div_col)))

MRD_div100_indel <- data.frame('div'=indel_col, 'MRD'=M100_MRDs_indel, 'read_len'=rep('100', length(indel_col)), 'type'=rep('Indels', length(indel_col)))
MRD_div150_indel <- data.frame('div'=indel_col, 'MRD'=M150_MRDs_indel, 'read_len'=rep('150', length(indel_col)), 'type'=rep('Indels', length(indel_col)))
MRD_div250_indel <- data.frame('div'=indel_col, 'MRD'=M250_MRDs_indel, 'read_len'=rep('250', length(indel_col)), 'type'=rep('Indels', length(indel_col)))
MRD_div500_indel <- data.frame('div'=indel_col, 'MRD'=M500_MRDs_indel, 'read_len'=rep('500', length(indel_col)), 'type'=rep('Indels', length(indel_col)))
MRD_div1000_indel <- data.frame('div'=indel_col, 'MRD'=M1000_MRDs_indel, 'read_len'=rep('1000', length(indel_col)), 'type'=rep('Indels', length(indel_col)))


MRD_div100_indel$MRD <- jitter(MRD_div100_indel$MRD, amount=2)
MRD_div150_indel$MRD <- jitter(MRD_div150_indel$MRD, amount=2)
MRD_div250_indel$MRD <- jitter(MRD_div250_indel$MRD, amount=2)
MRD_div500_indel$MRD <- jitter(MRD_div500_indel$MRD, amount=2)
MRD_div1000_indel$MRD <- jitter(MRD_div1000_indel$MRD, amount=2)

MRD_div <- rbind(MRD_div100, MRD_div150, MRD_div250, MRD_div500, MRD_div1000,
                 MRD_div100_indel, MRD_div150_indel, MRD_div250_indel, MRD_div500_indel, MRD_div1000_indel)

MRD_div_se <- summarySE(MRD_div, measurevar = 'MRD', groupvars = c('div', 'read_len', 'type'))

plot.mrd <- ggplot(MRD_div_se, aes(x=div, y=MRD, color=read_len)) + 
  facet_grid(. ~ type, scales='free') +
  #geom_point(data=MRD_div, alpha=0.3) +
  geom_point() +
  geom_errorbar(aes(ymin=MRD-se, ymax=MRD+se), width=0.7) +
  labs(color='Read length') + 
  ylab('error in coverage (%)') +
  xlab('simulated divergence (%)') +
  ggtitle('') +
  scale_color_manual(values=c('darkblue', 'orange', 'darkgoldenrod', 'darkgreen', 'darkred')) +
  geom_smooth(method='loess', span=.4, alpha=0.3, size=0.5, se=FALSE) +
  theme_minimal() 
#plot.mrd


# second plot about the correlation of observed vs expected divergence
DIV100 <- snp_ratio(rl100_list)
DIV150 <- snp_ratio(rl150_list)
DIV250 <- snp_ratio(rl250_list)
DIV500 <- snp_ratio(rl500_list)
DIV1000 <- snp_ratio(rl1000_list)

CORR_DIV100 <- data.frame('div'=div_col, 'perc_refsnp'=DIV100, 'read_len'=rep('100', length(div_col)))
CORR_DIV150 <- data.frame('div'=div_col, 'perc_refsnp'=DIV150, 'read_len'=rep('150', length(div_col)))
CORR_DIV250 <- data.frame('div'=div_col, 'perc_refsnp'=DIV250, 'read_len'=rep('250', length(div_col)))
CORR_DIV500 <- data.frame('div'=div_col, 'perc_refsnp'=DIV500, 'read_len'=rep('500', length(div_col)))
CORR_DIV1000 <- data.frame('div'=div_col, 'perc_refsnp'=DIV1000, 'read_len'=rep('1000', length(div_col)))

CORR_DIV <- rbind(CORR_DIV100, CORR_DIV150, CORR_DIV250, CORR_DIV500, CORR_DIV1000)
CORR_DIV_SE <- summarySE(CORR_DIV, measurevar = 'perc_refsnp', groupvars = c('div', 'read_len'))


plot.div <- ggplot(CORR_DIV_SE, aes(x=div, y=perc_refsnp, color=read_len)) + 
  #geom_point(data=CORR_DIV, alpha=0.3) +
  geom_point() +
  geom_errorbar(aes(ymin=perc_refsnp-ci, ymax=perc_refsnp+ci), width=.5) +
  geom_smooth(method='loess', span=.4, alpha=0.3, size=0.5, se=FALSE) +
  labs(color='Read length') + 
  ylab('observed (%)') +
  xlab('simulated (%)') +
  ggtitle('Divergence') +
  scale_color_manual(values=c('darkblue', 'orange', 'darkgoldenrod', 'darkgreen', 'darkred')) +
  scale_x_continuous(limits = c(0,27)) +
  theme_minimal() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=11))
#plot.div


