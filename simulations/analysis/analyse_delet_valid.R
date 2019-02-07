library(data.table)
library(gtools)
library(ggplot2)
library(summarySE)
library(cowplot)
# analyse internal deletion validation

sorted_file_list <- function(path, pattern){
  f <- list.files(path=path, pattern=pattern, full.names=TRUE, recursive=TRUE)
  f2 <- gsub(".", " ", f, fixed=TRUE)
  f2_sort <- mixedsort(f2)
  f2_sort_rep <- gsub(" ", ".", f2_sort, fixed=TRUE)
  return(f2_sort_rep)
}

extract_info <- function(file, rl){
  f_spl <- gsub("[_/.]", " ", file)
  f_spl2 <- trimws(gsub(pattern='[A-Za-z]', "", f_spl, perl=TRUE))
  f_l <- as.numeric(strsplit(f_spl2, split=' ')[[1]])
  f_vec <- f_l[!sapply(f_l, is.na)][-(1:1)]
  len <- f_vec[4] - f_vec[3]
    
  f <- fread(file)
  f_int_del <- f[!is.na(f$int_del), c(3,9,14)]
  
  if (nrow(f_int_del) != 0){
    int_del <- as.numeric(strsplit(f_int_del$int_del, ":")[[1]])
    start <- int_del[1]
    stop <- int_del[2]
    obs_cov <- int_del[3]
    mean_cov <- mean(f[(start + 1) :stop, cov])
  } else {
    obs_cov <- 0
    mean_cov <- 1
  }
  
  obs <- obs_cov / (obs_cov + mean_cov)
  sim <- f_vec[2] / 100
  file_res <- c(rl, sim, obs, len, obs_cov)
  return(file_res)
}

process_data <- function(path, rl){
  file_list <- sorted_file_list(path=path, pattern=".DOC5$") 
  M <- t(sapply(file_list, extract_info, rl))
  frame <- data.frame(M)
  names(frame) <- c('rl', 'sim', 'obs', 'len', 'obs_gap_cov')
  return(frame)
}

min_resid <- function(data, par){
  # find parameter combo to minimize residual sum of squares
  # intercept and slope fixed to prediction
  with(data, sum((0 + 1 * data[,"sim"] - (data[,"obs"]^par) )^2))
}

corr_fac <- function(frame){
  res <- optim(par=0, fn=min_resid, data=frame, method="Brent", lower=0, upper=5)
  return(res$par)
}

path <- '~/path/to/internal/deletion/sim/'

int_del80 <- process_data(path = paste(path, 'rl_80/', sep=''), rl = 80)
int_del100 <- process_data(path = paste(path, 'rl_100/', sep=''), rl = 100)
int_del175 <- process_data(path = paste(path, 'rl_175/', sep=''), rl = 175)
int_del250 <- process_data(path = paste(path, 'rl_250/', sep=''), rl = 250)
int_del500 <- process_data(path = paste(path, 'rl_500/', sep=''), rl = 500)
int_del1000 <- process_data(path = paste(path, 'rl_1000/', sep=''), rl = 1000)


int_del80$obs_corr <- int_del80$obs ^ corr_fac(frame=int_del80)
int_del100$obs_corr <- int_del100$obs ^ corr_fac(frame=int_del100)
int_del175$obs_corr <- int_del175$obs ^ corr_fac(frame=int_del175)
int_del250$obs_corr <- int_del250$obs ^ corr_fac(frame=int_del250)
int_del500$obs_corr <- int_del500$obs ^ corr_fac(frame=int_del500)
int_del1000$obs_corr <- int_del1000$obs ^ corr_fac(frame=int_del1000)

int_del <- rbind(int_del100, int_del250, int_del500, int_del1000)
int_del$rl <- as.factor(int_del$rl)
int_del$len <- as.factor(int_del$len)

int_del_se <- summarySE(int_del, measurevar = 'obs', groupvars = c('rl', 'sim'))
int_del_se_corr <- summarySE(int_del, measurevar = 'obs_corr', groupvars = c('rl', 'sim'))

# create dataframe for plotting points
raw <- int_del[,c(1,2,4,3)]
corr <- int_del[,c(1,2,4,6)]

plot.raw <- ggplot(int_del_se, aes(x=sim, y=obs, color=rl)) +
  geom_line(size=0.2) +
  geom_point() +
  geom_errorbar(aes(ymin=obs-se, ymax=obs+se), width=.01) +
  geom_jitter(data=raw, aes(shape=len), alpha=0.3, width=0.015, height=0.015) +
  labs(color='Read length', shape='Deletion size') + 
  ylab('observed freq.') +
  xlab('simulated freq.') +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  ggtitle('Raw') +
  scale_color_manual(values=c('darkblue', 'darkgoldenrod', 'darkgreen', 'darkred')) +
  scale_shape_manual(values=c(0,2,3,4,5,6,7)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=11))
#plot.raw

plot.corr <- ggplot(int_del_se_corr, aes(x=sim, y=obs_corr, color=rl)) +
  geom_line(size=0.2) +
  geom_point() +
  geom_jitter(data=corr, aes(shape=len), alpha=0.3, width=0.015, height=0.015) +
  labs(color='Read length', shape='Deletion size') + 
  ylab('estimated freq.') +
  xlab('simulated freq.') +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  ggtitle('Corrected') +
  scale_color_manual(values=c('darkblue', 'darkgoldenrod', 'darkgreen', 'darkred')) +
  scale_shape_manual(values=c(0,2,3,4,5,6,7)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=11))
#plot.corr

plot.leg <- get_legend(plot.raw)
plot.raw <- plot.raw + theme(legend.position="none")
plot.corr <- plot.corr + theme(legend.position="none")


# model the relation of read length and correction factor
corr_facs <- c(corr_fac(frame=int_del80),
               corr_fac(frame=int_del100),
               corr_fac(frame=int_del175),
               corr_fac(frame=int_del250),
               corr_fac(frame=int_del500),
               corr_fac(frame=int_del1000))

corr_facs <- as.data.frame(cbind(corr_facs, rl=c(80,100,175,250,500,1000)))

model <- lm(corr_facs$corr_facs ~ poly(corr_facs$rl,3, raw=TRUE))
#summary(model)
coefs <- unname(coef(model))
f <- function(x){  coefs[4]*x^3 + coefs[3]*x^2 + coefs[2]*x + coefs[1]   }

plot.model <- ggplot(corr_facs, aes(x=rl, y=corr_facs)) +
  stat_function(fun = f) +
  ylab("Corr. factor") +
  xlab("Read length") +
  scale_y_continuous(breaks=c(0.7,1.0,1.3)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1),
        text = element_text(size=8))






