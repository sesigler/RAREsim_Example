###############
#### Data was created/organized with Chromosome19_simulation_data.R
#### Create the plot with the final results
##### All of chromosome 19
#################

#### load in all the data:
##### Organized in the Chromosome19_simulation_data.R script
library(reshape2)
library(ggplot2)
library(gridExtra) # for grid.arrange
library(grid) # for textGrob and gpar
library(data.table) # for fwrite

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

Pop1 = 'AFR'
Pop2 = 'NFE'
Nsim = 20000
p_case = 120
p_conf_99 = 99
p_conf_95 = 95
p_conf_90 = 90
p_conf_80 = 80
bl <- 37  #block with median number of bases

dir_in <- 'C:/Users/sagee/Documents/GitHub/masters_project/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/'
# dir_out = 'C:/Users/sagee/Documents/GitHub/RAREsim_Example/'

# data frames to store all observed and expected values
obs_mac = c()
exp_mac = data.frame(matrix(nrow = 7, ncol = 9))

# Load in observed data (RAREsim)
fun_pcase = read.table(paste0(dir_in, 'obs_MACbin_fun', p_case, '.csv'), header = T, sep = ',')
fun_exp = read.table(paste0(dir_in, 'obs_MACbin_fun100.csv'), header = T, sep = ',')
fun_pconf_99 = read.table(paste0(dir_in, 'obs_MACbin_fun', p_conf_99, '.csv'), header = T, sep = ',')
fun_pconf_95 = read.table(paste0(dir_in, 'obs_MACbin_fun', p_conf_95, '.csv'), header = T, sep = ',')
fun_pconf_90 = read.table(paste0(dir_in, 'obs_MACbin_fun', p_conf_90, '.csv'), header = T, sep = ',')
fun_pconf_80 = read.table(paste0(dir_in, 'obs_MACbin_fun', p_conf_80, '.csv'), header = T, sep = ',')
syn_exp = read.table(paste0(dir_in, 'obs_MACbin_syn100.csv'), header = T, sep = ',')
syn_pconf_99 = read.table(paste0(dir_in, 'obs_MACbin_syn', p_conf_99, '.csv'), header = T, sep = ',')
syn_pconf_95 = read.table(paste0(dir_in, 'obs_MACbin_syn', p_conf_95, '.csv'), header = T, sep = ',')
syn_pconf_90 = read.table(paste0(dir_in, 'obs_MACbin_syn', p_conf_90, '.csv'), header = T, sep = ',')
syn_pconf_80 = read.table(paste0(dir_in, 'obs_MACbin_syn', p_conf_80, '.csv'), header = T, sep = ',')

# Load in expected data (gnomAD)
exp_fun_case = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_', p_case,  '.txt'), header=T, sep='\t')
exp_fun = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_100.txt'), header=T, sep='\t')
exp_fun_conf_99 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_', p_conf_99, '.txt'), header=T, sep='\t')
exp_fun_conf_95 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_', p_conf_95, '.txt'), header=T, sep='\t')
exp_fun_conf_90 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_', p_conf_90, '.txt'), header=T, sep='\t')
exp_fun_conf_80 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_fun_', p_conf_80, '.txt'), header=T, sep='\t')
exp_syn = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_syn_100.txt'), header=T, sep='\t')
exp_syn_conf_99 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_syn_', p_conf_99, '.txt'), header=T, sep='\t')
exp_syn_conf_95 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_syn_', p_conf_95, '.txt'), header=T, sep='\t')
exp_syn_conf_90 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_syn_', p_conf_90, '.txt'), header=T, sep='\t')
exp_syn_conf_80 = read.table(paste0(dir_in, 'MAC_bin_estimates_', Nsim, "_", Pop2, '_syn_', p_conf_80, '.txt'), header=T, sep='\t')

### CREATE TABLE of avg num of varints per bin for each pruning step
avg_vars <- data.frame(matrix(nrow = 11, ncol = 7))
colnames(avg_vars) <- c('Singletons', 'Doubletons', 'MAC=3-5', 'MAC=6-10', 'MAC=11-20',
                       'MAC=21-MAF=0.5%', 'MAF=0.5%-1%')
rownames(avg_vars) <- c('RAREsim functional-120%', 'RAREsim functional-100%', 'RAREsim functional-99%',
                        'RAREsim functional-95%', 'RAREsim functional-90%', 'RAREsim functional-80%',
                        'RAREsim synonymous-100%', 'RAREsim synonymous-99%', 'RAREsim synonymous-95%', 
                        'RAREsim synonymous-90%', 'RAREsim synonymous-80%')
avg_vars[1, ] <- lapply(fun_pcase[, 1:7], mean)
avg_vars[2, ] <- lapply(fun_exp[, 1:7], mean)
avg_vars[3, ] <- lapply(fun_pconf_99[, 1:7], mean)
avg_vars[4, ] <- lapply(fun_pconf_95[, 1:7], mean)
avg_vars[5, ] <- lapply(fun_pconf_90[, 1:7], mean)
avg_vars[6, ] <- lapply(fun_pconf_80[, 1:7], mean)
avg_vars[7, ] <- lapply(syn_exp[, 1:7], mean)
avg_vars[8, ] <- lapply(syn_pconf_99[, 1:7], mean)
avg_vars[9, ] <- lapply(syn_pconf_95[, 1:7], mean)
avg_vars[10, ] <- lapply(syn_pconf_90[, 1:7], mean)
avg_vars[11, ] <- lapply(syn_pconf_80[, 1:7], mean)

avg_vars_tpose <- data.frame(t(avg_vars))
colnames(avg_vars_tpose) <- rownames(avg_vars)

fwrite(avg_vars, paste0(dir_out, 'avg_variants_per_MACbin.csv'), quote=F, row.names=T, col.names=T, sep=',')
fwrite(avg_vars_tpose, paste0(dir_out, 'avg_variants_per_MACbin_tpose.csv'), quote=F, row.names=T, col.names=T, sep=',')

exp_all <- cbind(exp_fun_case, exp_fun[, 3], exp_fun_conf_99[, 3], exp_fun_conf_95[, 3],
                 exp_fun_conf_90[, 3], exp_fun_conf_80[, 3], exp_syn[, 3], exp_syn_conf_99[, 3], 
                 exp_syn_conf_95[, 3], exp_syn_conf_90[, 3], exp_syn_conf_80[, 3])

colnames(exp_all) <- c('Lower', 'Upper', 'Expected_var-Fun 120%', 'Expected_var-Fun 100%',
                       'Expected_var-Fun 99%', 'Expected_var-Fun 95%', 'Expected_var-Fun 90%',
                       'Expected_var-Fun 80%','Expected_var-Syn 100%', 'Expected_var-Syn 99%',
                       'Expected_var-Syn 95%', 'Expected_var-Syn 90%', 'Expected_var-Syn 80%')

fwrite(exp_all, paste0(dir_out, 'exp_variants_per_MACbin_20000_NFE_all.csv'), quote=F, row.names=F, col.names=T, sep=',')
##################

# Form the data frames for the observed and expected counts
obs_mac <- rbind(fun_pcase, fun_exp, fun_pconf1, fun_pconf2, syn_exp, syn_pconf1, syn_pconf2)
colnames(obs_mac) <- c('Singletons', 'Doubletons', 'MAC=3-5', 'MAC=6-10', 'MAC=11-20',
                       'MAC=21-MAF=0.5%', 'MAF=0.5%-1%', 'rep', 'data')

exp_mac[1:7, 1:7] <- rbind(exp_fun_case$Expected_var, exp_fun$Expected_var, 
                           exp_fun_conf1$Expected_var, exp_fun_conf2$Expected_var,
                           exp_syn$Expected_var, exp_syn_conf1$Expected_var, 
                           exp_syn_conf2$Expected_var)
colnames(exp_mac) <- colnames(obs_mac)
exp_mac$rep <- '.'
exp_mac$data <- c('gnomAD functional-120%', 'gnomAD functional-100%', 
                  paste0('gnomAD functional-', p_conf1, '%'), 
                  paste0('gnomAD functional-', p_conf2, '%'),
                  'gnomAD synonymous-100%', paste0('gnomAD synonymous-', p_conf1, '%'),
                  paste0('gnomAD synonymous-', p_conf2, '%'))

# Combine the observed and expected data into one dataframe
melted_obs <- reshape2::melt(obs_mac, id = c('data', 'rep'))
melted_exp <- reshape2::melt(exp_mac, id = c('data', 'rep'))

# setwd('/Users/megansorenson/Documents/Package/RAREsim_Example/Code_from_analysis/Simulation/')
# setwd('C:/Users/sagee/Documents/GitHub/RAREsim_Example/Code_from_analysis/Simulation/')
# mac <- read.table('Simulation_results100reps_chr19_final.txt', header = TRUE, sep = '\t')
# mac_all <- read.table('mac_all_chr19_final.txt', header = TRUE, sep = '\t')
# melted_mac <- melt(mac, id = c('block', 'pop', 'data', 'rep'))
# melted_mac_all <- melt(mac_all, id = c('block', 'pop', 'data', 'rep'))
# head(melted_mac)
# table(melted_mac$pop, melted_mac$data)


### Separate gnomAD and the simulation data
# melted_gnom <- melted_mac[which(melted_mac$data == 'gnomAD'),]
# melted_mac1 <- melted_mac[-c(which(melted_mac$data == 'gnomAD')),]

# melted_gnom_all <- melted_mac_all[which(melted_mac_all$data == 'gnomAD functional-120%' | melted_mac_all$data == 'gnomAD functional-100%' | melted_mac_all$data == 'gnomAD functional-99%'),]
# melted_mac1_all <- melted_mac_all[-c(which(melted_mac_all$data == 'gnomAD' | melted_mac_all$data == 'gnomAD functional' | melted_mac_all$data == 'gnomAD synonymous')),]

# melted_mac1$data <- as.character(melted_mac1$data)
# melted_mac1$data[which(melted_mac1$data  == 'HAPGEN2 over-simulated')] <- 'HAPGEN2 with all bp'
# melted_mac1$data[which(melted_mac1$data  == 'Original HAPGEN2')] <- 'HAPGEN2 with polymorphic SNVs'
# 
# melted_mac1_all$data <- as.character(melted_mac1_all$data)
# melted_mac1_all$data[which(melted_mac1_all$data  == 'HAPGEN2 over-simulated')] <- 'HAPGEN2 with all bp'
# melted_mac1_all$data[which(melted_mac1_all$data  == 'Original HAPGEN2')] <- 'HAPGEN2 with polymorphic SNVs'


cbPalette <- c("#56B4E9","#E69F00",  
               "#009E73", "#CC79A7", "#0072B2",
               "#D55E00")

cbPalette <- c("#56B4E9", "#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7")

cols_fun <- c("#009E73", "#CC79A7") #Functional
cols_syn <- c("#0072B2", "#D55E00") #Synonymous
cols_all <- c("#009E73", "#0072B2", "#CC79A7", "#D55E00")

# Functional-120%
fun120 <- ggplot(subset(melted_obs, data %in% paste0('RAREsim functional-', p_case, '%')), 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = paste0('Observed vs Expected Number of Variants per MAC bin \nPop: 100% NFE \nPruning: ', p_case, '% Functional'), 
       y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% paste0('gnomAD functional-', p_case, '%')),
             aes(x=variable, y=value), shape  = 8) +
  # theme(axis.title.x = element_blank()) +
  # theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_fun) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.2, 0.5, 0.1, 0.2), "cm"))
fun120
ggsave(file = paste0(dir_out, 'obs_v_exp_MACbins_', Pop2, '_fun', p_case, '.jpg'),
       plot = fun120, height = 8, width = 12, units = 'in')

# Functional-100%
fun100 <- ggplot(subset(melted_obs, data %in% 'RAREsim functional-100%'), 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = 'Pruning: 100% Functional', y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% 'gnomAD functional-100%'),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_fun) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
fun100

# Synonymous-100%
syn100 <- ggplot(subset(melted_obs, data %in% 'RAREsim synonymous-100%'), 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = 'Pruning: 100% Synonymous', y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% 'gnomAD synonymous-100%'),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_syn) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
syn100

# Functional-p_conf1 %
fun_p_conf1 <- ggplot(subset(melted_obs, data %in% paste0('RAREsim functional-', p_conf1, '%')), 
                 aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = paste0('Pruning: ', p_conf1, '% Functional'), 
       y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% paste0('gnomAD functional-', p_conf1, '%')),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_fun) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
fun_p_conf1

# Synonymous-p-conf1 %
syn_p_conf1 <- ggplot(subset(melted_obs, data %in% paste0('RAREsim synonymous-', p_conf1, '%')), 
                 aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = paste0('Pruning: ', p_conf1, '% Synonymous'), 
       y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% paste0('gnomAD synonymous-', p_conf1, '%')),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_syn) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
syn_p_conf1

# Functional-p_conf2 %
fun_p_conf2 <- ggplot(subset(melted_obs, data %in% paste0('RAREsim functional-', p_conf2, '%')), 
                aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = paste0('Pruning: ', p_conf2, '% Functional'), 
       y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% paste0('gnomAD functional-', p_conf2, '%')),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_fun) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
fun_p_conf2

# Synonymous-p_conf2 %
syn_p_conf2 <- ggplot(subset(melted_obs, data %in% paste0('RAREsim synonymous-', p_conf2, '%')), 
                aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(title = paste0('Pruning: ', p_conf2, '% Synonymous'), 
       y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% paste0('gnomAD synonymous-', p_conf2, '%')),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_syn) +
  # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
syn_p_conf2

### SAVE THE PLOTS ###
# Extract the legend for each pruning scenario
both100 <- ggplot(subset(melted_obs, data %in% c('RAREsim functional-100%', 'RAREsim synonymous-100%')), 
                 aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% c('gnomAD functional-100%', 'gnomAD synonymous-100%')),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_all) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))

both_p_conf1 <- ggplot(subset(melted_obs, data %in% c(paste0('RAREsim functional-', p_conf1, '%'), 
                                                      paste0('RAREsim synonymous-', p_conf1, '%'))), 
                  aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% c(paste0('gnomAD functional-', p_conf1, '%'), 
                                                 paste0('gnomAD synonymous-', p_conf1, '%'))),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_all) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))

both_p_conf2 <- ggplot(subset(melted_obs, data %in% c(paste0('RAREsim functional-', p_conf2, '%'),
                                                      paste0('RAREsim synonymous-', p_conf2, '%'))), 
                 aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bins') +
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=subset(melted_exp, data %in% c(paste0('gnomAD functional-', p_conf2, '%'),
                                                 paste0('gnomAD synonymous-', p_conf2, '%'))),
             aes(x=variable, y=value), shape  = 8) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=cols_all) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))

leg100 <- g_legend(both100)
leg_p_conf1 <- g_legend(both_p_conf1)
leg_p_conf2 <- g_legend(both_p_conf2)

# Combine the Functional and Synonymous plots into one figure
pruned100 <- grid.arrange(arrangeGrob(fun100 + theme(legend.position="none"),
                                      syn100 + theme(legend.position="none"),
                                      nrow=1),
                         leg100, nrow=2, heights = c(10, 1),
                         left = 'Number of Variants',
                         top = textGrob('100% NFE: Observed vs Expected Number of Variants per MAC bin',
                                        gp = gpar(fontsize = 18, font = 1)))
# Save the plot
ggsave(file = paste0(dir_out, 'obs_v_exp_MACbins_', Pop2, '_fun100_syn100.jpg'),
       plot = pruned100, height = 8, width = 12, units = 'in')

pruned_p_conf1 <- grid.arrange(arrangeGrob(fun_p_conf1 + theme(legend.position="none"),
                                           syn_p_conf1 + theme(legend.position="none"),
                                           nrow=1),
                          leg_p_conf1, nrow=2, heights = c(10, 1),
                          left = 'Number of Variants',
                          top = textGrob('100% NFE: Observed vs Expected Number of Variants per MAC bin',
                                         gp = gpar(fontsize = 18, font = 1)))
# Save the plot
ggsave(file = paste0(dir_out, 'obs_v_exp_MACbins_', Pop2, '_fun', p_conf1, '_syn', p_conf1, '.jpg'),
       plot = pruned_p_conf1, height = 8, width = 12, units = 'in')

pruned_p_conf2 <- grid.arrange(arrangeGrob(fun_p_conf2 + theme(legend.position="none"),
                                           syn_p_conf2 + theme(legend.position="none"),
                                           nrow=1),
                         leg_p_conf2, nrow=2, heights = c(10, 1),
                         left = 'Number of Variants',
                         top = textGrob('100% NFE: Observed vs Expected Number of Variants per MAC bin',
                                        gp = gpar(fontsize = 18, font = 1)))
# Save the plot
ggsave(file = paste0(dir_out, 'obs_v_exp_MACbins_', Pop2, '_fun', p_conf2, '_syn', p_conf2, '.jpg'),
       plot = pruned_p_conf2, height = 8, width = 12, units = 'in')


# for(bl in c(56,66)){ # 5th and 95th percentile for number of bases
#   mm <- melted_mac1[which(melted_mac1$block  == bl),]
#   gm <- melted_gnom[which(melted_gnom$block  == bl),]
#   
#   
#   pop <- 'AFR'
#   p1 <- ggplot(mm[which(mm$pop == pop),], 
#                aes(x=variable, y=value , col = data)) +
#     geom_boxplot() + 
#     labs(y = 'Number of Variants', x = 'MAC Bin')+
#     theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
#     geom_point(data=gm[which(gm$pop == pop),],
#                aes(x=variable, y=value), shape  = 8)+
#     theme(legend.position="bottom" )+
#     theme(legend.title = element_blank()) +
#     theme(axis.title.x = element_blank())+
#     scale_color_manual(values=cbPalette)+
#     theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
#   p1
#   
#   pop <- 'EAS'
#   p2 <- ggplot(mm[which(mm$pop == pop),], 
#                aes(x=variable, y=value , col = data)) +
#     geom_boxplot() + 
#     labs(y = 'Number of Variants', x = 'MAC Bin')+
#     theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
#     geom_point(data=gm[which(gm$pop == pop),],
#                aes(x=variable, y=value), shape  = 8)+
#     theme(axis.title.x = element_blank())+
#     theme(axis.title.y = element_blank())+
#     scale_color_manual(values=cbPalette)+
#     theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
#   p2
#   
#   pop <- 'NFE'
#   p3 <- ggplot(mm[which(mm$pop == pop),], 
#                aes(x=variable, y=value , col = data)) +
#     geom_boxplot() + 
#     labs(y = 'Number of Variants', x = 'MAC Bin')+
#     theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
#     geom_point(data=gm[which(gm$pop == pop),],
#                aes(x=variable, y=value), shape  = 8)+
#     theme(axis.title.x = element_blank())+
#     theme(axis.title.y = element_blank())+
#     scale_color_manual(values=cbPalette)+
#     theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
#   p3
#   
#   pop <- 'SAS'
#   p4 <- ggplot(mm[which(mm$pop == pop),], 
#                aes(x=variable, y=value , col = data)) +
#     geom_boxplot() + 
#     labs(y = 'Number of Variants', x = 'MAC Bin')+
#     theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
#     geom_point(data=gm[which(gm$pop == pop),],
#                aes(x=variable, y=value), shape  = 8)+
#     theme(axis.title.x = element_blank())+
#     theme(axis.title.y = element_blank())+
#     scale_color_manual(values=cbPalette)+
#     theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
#   p4 
#   
#   
#   #####  try to extract the legend:
#   g_legend<-function(a.gplot){
#     tmp <- ggplot_gtable(ggplot_build(a.gplot))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     return(legend)}
#   
#   mylegend<-g_legend(p1)
#   
#   p5 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
#                                  p2 + theme(legend.position="none"),
#                                  p3 + theme(legend.position="none"),
#                                  p4 + theme(legend.position="none"),
#                                  nrow=1),
#                      mylegend, nrow=2,heights=c(10, 1.1))
#   
#   # ggsave(file = paste0('Block',bl,'main.jpg'),
#   #        plot = p5, height = 5, width = 10, units = 'in')
# }


for(bl in c(37)){ ## The median block
  # mm <- melted_mac1[which(melted_mac1$block  == bl),]
  # gm <- melted_gnom[which(melted_gnom$block  == bl),]
  # mm <- melted_mac1_all[which(melted_mac1_all$block  == bl),]
  # gm <- melted_gnom_all[which(melted_gnom_all$block  == bl),]
  
  pop <- 'AFR'
  p1 <- ggplot(mm[which(mm$pop == pop),], 
               aes(x=variable, y=value , col = data)) +
    geom_boxplot() + 
    labs(y = 'Number of Variants', x = 'MAC Bin')+
    theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
    geom_point(data=gm[which(gm$pop == pop),],
               aes(x=variable, y=value), shape  = rep(c(8, 1, 4), 7)) +
    theme(legend.position="bottom")+
    theme(legend.title = element_blank()) +
    theme(axis.title.x = element_blank()) +
    scale_color_manual(values=cbPalette) +
    guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
    theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
  p1
  
  # pop <- 'EAS'
  # p2 <- ggplot(mm[which(mm$pop == pop),], 
  #              aes(x=variable, y=value , col = data)) +
  #   geom_boxplot() + 
  #   labs(y = 'Number of Variants', x = 'MAC Bin')+
  #   theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  #   geom_point(data=gm[which(gm$pop == pop),],
  #              aes(x=variable, y=value), shape  = 8)+
  #   theme(axis.title.x = element_blank())+
  #   scale_color_manual(values=cbPalette)+
  #   theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
  # p2
  
  # pop <- 'NFE'
  p3 <- ggplot(melted_obs[which(melted_obs$data == 'RAREsim functional-120%'),], 
               aes(x=variable, y=value , col = data)) +
    geom_boxplot() + 
    labs(y = 'Number of Variants', x = 'MAC Bin') +
    theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
    geom_point(data=melted_exp[which(melted_exp$data == 'gnomAD functional-120%'),],
               aes(x=variable, y=value), shape  = 8) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    # scale_color_manual(values=cbPalette) +
    # guides(color = guide_legend(override.aes=list(shape = c(8, 1, 4, 8, 8, 8)))) +
    theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
  p3
  
  # pop <- 'SAS'
  # p4 <- ggplot(mm[which(mm$pop == pop),], 
  #              aes(x=variable, y=value , col = data)) +
  #   geom_boxplot() + 
  #   labs(y = 'Number of Variants', x = 'MAC Bin')+
  #   theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  #   geom_point(data=gm[which(gm$pop == pop),],
  #              aes(x=variable, y=value), shape  = 8)+
  #   theme(axis.title.x = element_blank())+
  #   theme(axis.title.y = element_blank())+
  #   scale_color_manual(values=cbPalette)+
  #   theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
  # p4 
  
  
  ##### Extract the legend:
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
    }
  
  mylegend<-g_legend(p1)
  
  p5 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                 p3 + theme(legend.position="none"),
                                 nrow=1),
                     mylegend, nrow=2,heights=c(10, 1.1))
  
  # ggsave(file = paste0(dir_out, 'Block',bl,'main_2only.jpg'),
  #        plot = p5, height = 5, width = 8, units = 'in')
  ggsave(file = paste0(dir_out, 'Block',bl,'main_2only_stratified.jpg'),
         plot = p5, height = 5, width = 8, units = 'in')
  
  p6 <- grid.arrange(arrangeGrob(p2 + theme(legend.position="none") ,
                                 p4 + theme(legend.position="none"),
                                 nrow=1),
                     mylegend, nrow=2,heights=c(10, 1.1))
  
  # ggsave(file = paste0('Block',bl,'supp_2only.jpg'),
  #        plot = p6, height = 5, width = 8, units = 'in')
  
}


###### Look at all of chromosome 19

dt <- levels(as.factor(mac$data))
head(mac)
all <- c()
for(pop in c('AFR', 'EAS',  'NFE', 'SAS')){
  df <- mac[which(mac$pop == pop),]
  for(i in 1:length(dt)){
    dff <- df[which(df$data == dt[i]),]
    for(rep in 1:100){
      df1 <- dff[which(dff$rep == rep),]
      temp <- colSums(df1[,3:9])
      temp <- as.data.frame(t(temp))
      temp$data <- dt[i]
      temp$pop <- pop
      temp$rep <- rep
      all  <- rbind(all, temp)
    }
    
  }
  dff <- df[which(df$data == 'gnomAD'),]
  temp <- colSums(dff[,3:9])
  temp <- as.data.frame(t(temp))
  temp$data <-  'gnomAD'
  temp$pop <- pop
  temp$rep <- '.'
  all  <- rbind(all,  temp)
}

all <- all[-c(which(all$data == 'gnomAD' & all$Singletons == 0 )),]


library(reshape2)
all_melt <- melt(all, id = c('pop', 'data', 'rep'))


all_melt1 <- all_melt[-c(which(all_melt$data == 'gnomAD')),]
gnom_melt <- all_melt[which(all_melt$data == 'gnomAD'),]

table(gnom_melt$pop, gnom_melt$variable)
table(all_melt1$data)

pop <- 'AFR'
p1 <- ggplot(all_melt1[which(all_melt1$pop == pop),], 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bin')+
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=gnom_melt[which(gnom_melt$pop == pop),],
             aes(x=variable, y=value), shape  = 8)+
  theme(legend.position="bottom" )+
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())+
  scale_color_manual(values=cbPalette)+
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
p1

pop <- 'EAS'
p2 <- ggplot(all_melt1[which(all_melt1$pop == pop),], 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bin')+
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=gnom_melt[which(gnom_melt$pop == pop),],
             aes(x=variable, y=value), shape  = 8)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values=cbPalette)+
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
p2

pop <- 'NFE'
p3 <- ggplot(all_melt1[which(all_melt1$pop == pop),], 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bin')+
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=gnom_melt[which(gnom_melt$pop == pop),],
             aes(x=variable, y=value), shape  = 8)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values=cbPalette)+
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
p3

pop <- 'SAS'
p4 <- ggplot(all_melt1[which(all_melt1$pop == pop),], 
             aes(x=variable, y=value , col = data)) +
  geom_boxplot() + 
  labs(y = 'Number of Variants', x = 'MAC Bin')+
  theme(axis.text.x = element_text(angle = 35, hjust=0.65)) +
  geom_point(data=gnom_melt[which(gnom_melt$pop == pop),],
             aes(x=variable, y=value), shape  = 8)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values=cbPalette)+
  theme(plot.margin = unit(c(0.2,0.5,0.1,0.2), "cm"))
p4 


#####  try to extract the legend:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

p5 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               p3 + theme(legend.position="none"),
                               p4 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1.1))
# 
# ggsave(file ='All_chr19_main.jpg',
#        plot = p5, height = 5, width = 10, units = 'in')





