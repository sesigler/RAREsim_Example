###############
#### Organize the simulation data
##### All of chromosome 19
##### To visualize, run Chromosome19_simulations_boxplots.R
#################

## gnomAD data:
# setwd('/Users/megansorenson/Documents/Package/RAREsim_Example/Code_from_analysis/Simulation/')
setwd('C:/Users/sagee/OneDrive/Documents/GitHub/RAREsim_Example/Stratified_Target_Data/')
gnom_strat <- read.table('AFS_target_data_stratified.txt', header =  TRUE, sep  = '\t')
nvar_gnom_strat <- read.table('Nvariants_target_data_stratified.txt', header =  TRUE, sep  = '\t')

library(data.table)
setwd('C:/Users/sagee/OneDrive/Documents/GitHub/RAREsim_Example/Reference_Data/')


### Convert AFs to ACs
# Note that the AFS is only for rare variants which is why the AFs don't sum to 1
# or the ACs don't sume to the total
gnom_strat$Singletons <- gnom_strat$MAC_1*gnom_strat$total
gnom_strat$Doubletons <- gnom_strat$MAC_2*gnom_strat$total
gnom_strat$MAC.3.5 <- gnom_strat$MAC_3.5*gnom_strat$total
gnom_strat$MAC.6.10 <- gnom_strat$MAC_6.10*gnom_strat$total
gnom_strat$MAC.11.20 <- gnom_strat$MAC_11.20*gnom_strat$total
gnom_strat$MAC.21.MAF005 <- gnom_strat$MAC_21.MAF0.5*gnom_strat$total
gnom_strat$MAF005.MAF01 <- gnom_strat$MAF0.5.1*gnom_strat$total


setwd('C:/Users/sagee/OneDrive/Documents/GitHub/RAREsim_Example/Code_from_analysis/Simulation/')
gnom <- read.table('AFS_data_blocks_merged.txt', header =  TRUE, sep  = '\t')


#### Expected variants based on the functions we fit
##### This is the goal that RAREsim will prune to:
fit <- read.table('Expected_counts4raresim_final.txt', header = TRUE)
head(fit)

block <- levels(as.factor(as.character(fit$block)))

###### Load in the original simulation data

# ## Create a template to load in each bin:
# mac_template <- as.data.frame(matrix(nrow=100, ncol = 10))
# colnames(mac_template) <- c('block', 'pop', 'Singletons', 'Doubletons', 'MAC=3-5',
#                             'MAC=6-10', 'MAC=11-20', 
#                             'MAC=21-MAF005',
#                             'MAF005-MAF01', 'rep')
# 
# #### loop defines the bins
# loop <- as.data.frame(matrix(nrow=7,  ncol = 2))
# colnames(loop) <- c('Start', 'End')
# loop$Start <- c(1,2,3,6,11,21,'.')
# loop$End <- c(1,2,5,10,20,'.','.')
### N will change with  each ancestry now
# 
#### mac1 is the original hapgen2
# mac1 <- mac_template
# 
# setwd('/Users/megansorenson/Documents/RAREsim/Chr19/Simulations_gnomADsize/OG/')
# mac1 <- c()
# for(pop in c('EAS', 'SAS', 'AFR', 'NFE')){ # 
#   N <- fit$N[which(fit$pop == pop)]
#   N005 <- floor(0.005*2*N)
#   loop$End[6] <- N005
#   loop$Start[7] <- N005 + 1
#   N01 <- floor(0.01*2*N)
#   loop$End[7] <- N01
#   loop$Start <- as.numeric(loop$Start)
#   loop$End <- as.numeric(loop$End)
#   for(bl in block){ #
#     mac <- mac_template
#     mac$block <- bl
#     mac$pop <- pop
#     for(rep in c(1:100)){
#       temp <- read.table(paste0('MAC_',pop,'/Block',bl,'/Rep', rep,'.count.txt'))
#       
#       for(k in 1:nrow(loop)){
#         tmp <-  temp[which(temp >= loop$Start[k] & temp <= loop$End[k]),]
#         n1 <- length(tmp)
#         mac[rep,k+2] <- n1
#       }
#       mac$rep[rep] <- rep
#     }
#     mac1 <- rbind(mac1, mac)
#   }
# }
# 
# mac1$data <- 'Original HAPGEN2'
# 
# write.table(mac1, 'All_original_hapgen2_AC_information.txt',
#             row.names  = FALSE, sep = '\t', quote = FALSE)

mac1 <-  read.table('All_original_hapgen2_AC_information.txt', header=TRUE,
                    sep = '\t')

##############
### load  in the over-simulated allele counts
################

# mac2 <- c()
# for(pop in c('AFR', 'EAS', 'SAS', 'NFE')){
#   mac <- mac_template
#   N <- fit$N[which(fit$pop == pop)]
#   N005 <- floor(0.005*2*N)
#   loop$End[6] <- N005
#   loop$Start[7] <- N005 + 1
#   N01 <- floor(0.01*2*N)
#   loop$End[7] <- N01
#   loop$Start <- as.numeric(loop$Start)
#   loop$End <- as.numeric(loop$End)
#   for(bl in block){
#     mac <- mac_template
#     mac$block <- bl
#     mac$pop <- pop
#     for(rep in c(1:100)){
#       temp <- read.table(paste0(pop,'_MAC/Block',bl,'/Rep', rep,'.count.txt'))
#       for(k in 1:nrow(loop)){
#         tmp <-  temp[which(temp >= loop$Start[k] & temp <= loop$End[k]),]
#         n1 <- length(tmp)
#         mac[rep,k+2] <- n1
#       }
#       mac$rep[rep] <- rep
#     }
# 
#     mac2 <- rbind(mac2, mac)
#   }
#   print(pop)
#   }
# 
# mac2$data <- 'HAPGEN2 over-simulated'
# 
# write.table(mac2, 'All_oversimulated_hapgen2_AC_information.txt',
#            row.names  = FALSE, sep = '\t', quote = FALSE)


mac2 <- read.table('All_oversimulated_hapgen2_AC_information.txt', 
                   header  =  TRUE, sep= '\t')

#############
##### RAREsim simulated data:
############

# mac3 <- c()
# set.seed(1234)
# for(pop in c('EAS', 'SAS', 'AFR', 'NFE')){
#   N <- fit$N[which(fit$pop == pop)]
#   N005 <- floor(0.005*2*N)
#   loop$End[6] <- N005
#   loop$Start[7] <- N005 + 1
#   N01 <- floor(0.01*2*N)
#   loop$End[7] <- N01
#   loop$Start <- as.numeric(loop$Start)
#   loop$End <- as.numeric(loop$End)
#   for(bl in  block){
#     mac <- mac_template
#     mac$block <- bl
#     mac$pop <- pop
# 
#     for(rep in c(1:100)){
#       rem_all<-c()
#       change_all <-  c()
#       temp <- read.table(paste0(pop,'_MAC/Block',bl,'/Rep', rep,'.count.txt'))
#       temp <- as.data.frame(temp)
#       temp$num <-  1:nrow(temp)
#       for(k in nrow(loop):1){
#         tmp <-  temp[which(temp$V1 >= loop$Start[k] & temp$V1 <= loop$End[k]),]
#         ##### if there are not any SNVs in the bin, tmp is all NA's...
#         n1 <- nrow(tmp)
#         expect <- fit_all[which(fit_all$pop == pop & fit_all$block == bl),k+2]  ##expected
#         #print(expect)
#         if(n1>0){
#           p <- expect/n1
#         }else{
#           if(expect < 1){next}
#           if(expect >=1){
#             p <- expect
#           }
#         }
#         if(p <= 1){
# 
#           rd <- runif(n1)
# 
#           rem <- tmp[c(which(rd >= p)),]
# 
#           #print(rem)
#           ### now need to keep these ones...
# 
#           rem_all <- rbind(rem_all, rem)
# 
#         }
#         if(p>1){
#           if(k==7){
#             if((expect-n1) >2){
#               print(c('Error!', pop, bl,n1, expect, (n1-expect)))}
#             next}
#           if(length(rem_all)>0){
#           n2 <- expect  - n1 ### n2 is how many more we need
#           p1 <- n2/nrow(rem_all)  ### this is the proportion from rem_all
#           rem_all$rd <- runif(nrow(rem_all))
#           to_change <- rem_all[which(rem_all$rd  <= p1),]
#           mac_range <-  c(loop$Start[k]:loop$End[k])
#           to_change$new_MAC <- sample(mac_range, nrow(to_change),
#                                       replace = TRUE)
#           ### change them within the MAC file.
#           ####  this needs to change for the haplotype file!!!
#           temp$V1[which(temp$num %in% to_change$num)]<- to_change$new_MAC
#           rem_all <- rem_all[which(rem_all$rd > p1),]
#           change_all<- rbind(change_all,  to_change)
#           }else{
#             print('Error2')
#             next}
#         }
#       }
#       temp <- temp[-c(rem_all$num),]
# 
#       ##### At this point, set rem_all to all reference alleles in haplotype files
#       ##### And set change_all to the given MAC
# 
#       #####  ok, now I need the allele count to add to the plot
#       for(k in 1:nrow(loop)){
#         tmp <-  temp[which(temp$V1 >= loop$Start[k] & temp$V1 <= loop$End[k]),]
#         n1 <- nrow(tmp)
#         mac[rep,k+2] <- n1
#       }
#       mac$rep[rep] <- rep
#     }
#     mac3  <- rbind(mac3,mac)
# 
#   }
# }
# 
# 
# mac3$data  <- 'RAREsim'
# 
# write.table(mac3, 'All_RAREsim_AC_information.txt',
#             row.names  = FALSE, sep = '\t', quote = FALSE)

mac3 <- read.table('All_RAREsim_AC_information_final.txt', header = TRUE,
                   sep = '\t')

### Bind them all for the graph/summary
names(mac1) <- names(mac2)
names(mac3) <- names(mac2)
mac <- rbind(mac1,mac2)
names(mac)
names(mac3)
mac <- rbind(mac, mac3)

#### add gnomAD to the plot summary
names(mac)
names(gnom)
gnom <- gnom[,c(1,10,2:8)]
gnom$rep <- '.'
gnom$data <- 'gnomAD'

####
names(gnom_strat)
gnom_strat <- gnom_strat[, c(1, 10, 12:18, 11)]
gnom_strat$rep <- '.'
gnom_fun <- gnom_strat[gnom_strat$fun == 'fun', ]
gnom_syn <- gnom_strat[gnom_strat$fun == 'syn', ]
gnom_fun$data <- 'gnomAD functional'
gnom_syn$data <- 'gnomAD synonymous'
gnom_fun <- subset(gnom_fun, select = -fun)
gnom_syn <- subset(gnom_syn, select = -fun)
names(gnom) <- names(mac)

mac_all <- rbind(mac, gnom)
mac_all <- rbind(mac_all, gnom_fun)
mac_all <- rbind(mac_all, gnom_syn)
####

names(gnom) <- names(mac)
mac <- rbind(mac,gnom)

head(mac)

### Pretty names
colnames(mac)[5] <- 'MAC=3-5'
colnames(mac)[6] <- 'MAC=6-10'
colnames(mac)[7] <- 'MAC=11-20'
colnames(mac)[8] <- 'MAC=21-MAF=0.5%'
colnames(mac)[9] <- 'MAF=0.5%-1%'

table(mac$data)
# write.table(mac, 'Simulation_results100reps_chr19_final.txt', row.names = FALSE,
#             quote = FALSE, sep = '\t')

colnames(mac_all)[5] <- 'MAC=3-5'
colnames(mac_all)[6] <- 'MAC=6-10'
colnames(mac_all)[7] <- 'MAC=11-20'
colnames(mac_all)[8] <- 'MAC=21-MAF=0.5%'
colnames(mac_all)[9] <- 'MAF=0.5%-1%'

table(mac_all$data)
write.table(mac_all, 'mac_all_chr19_final.txt', row.names = FALSE,
            quote = FALSE, sep = '\t')


##### get the results for the final table:

for(pop in c('AFR', 'EAS', 'NFE','SAS')){

temp <- melted_mac[which(melted_mac$pop  == pop),]
keep <- as.data.frame(matrix(nrow = 707, ncol = 12 ))
colnames(keep) <- c('Block', 'MAC Bin', 'RAREsim Median', 'RAREsim Minimum', 'RAREsim Maximum',
                    'Original HAPGEN2 Median', 'Original HAPGEN2 Minimim', 'Original HAPGEN2 Maximum',
                     'HAPGEN2 with All bp Median', 'HAPGEN2 with All bp  Minimim', 'HAPGEN2 with All bp  Maximum')

keep[,c(1,2)] <- temp[which(temp$data == 'RAREsim' & temp$rep == 1),c(1,5)]

colnames(keep) <- c('Block', 'MAC Bin', 'gnomAD', 'RAREsim Median', 'RAREsim Minimum', 'RAREsim Maximum',
                    'Original HAPGEN2 Median', 'Original HAPGEN2 Minimim', 'Original HAPGEN2 Maximum',
                    'HAPGEN2 with All bp Median', 'HAPGEN2 with All bp  Minimim', 'HAPGEN2 with All bp  Maximum')

for(i in 1:nrow(keep)){
  tmp <- temp[which(temp$block == keep$Block[i] & temp$variable == keep$`MAC Bin`[i]),]
  keep[i,3] <- median(tmp$value[which(tmp$data == 'gnomAD')])

  keep[i,4] <- median(tmp$value[which(tmp$data == 'RAREsim')])
  keep[i,5] <- min(tmp$value[which(tmp$data == 'RAREsim')])
  keep[i,6] <- max(tmp$value[which(tmp$data == 'RAREsim')])

  keep[i,7] <- median(tmp$value[which(tmp$data == 'Original HAPGEN2')])
  keep[i,8] <- min(tmp$value[which(tmp$data == 'Original HAPGEN2')])
  keep[i,9] <- max(tmp$value[which(tmp$data == 'Original HAPGEN2')])

  keep[i,10] <- median(tmp$value[which(tmp$data == 'HAPGEN2 over-simulated')])
  keep[i,11] <- min(tmp$value[which(tmp$data == 'HAPGEN2 over-simulated')])
  keep[i,12] <- max(tmp$value[which(tmp$data == 'HAPGEN2 over-simulated')])
}

# write.table(keep, paste0(pop, '_simulation_result_summary.txt'),
#                          sep = '\t', quote = FALSE,
#             row.names = FALSE)
}

