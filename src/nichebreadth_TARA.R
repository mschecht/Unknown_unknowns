################################################
# Levin's Niche Breadth Calculation for Clusters
################################################

##############################################################
# R data sets that calculated in this script:
# - component abundance data filtered for 1e-5 mean proportion
##############################################################

# To reproduce this workflow upload clstrs_tara_f2.RData from Figshare (10.6084/m9.figshare.5979658).

#################
# Set environment
#################
library(tidyverse)
library(pbmcapply)
library(spaa)
library(vegan)
library(rworldmap)
library(RSQLite)
library(ggpubr)
library(data.table)

######################################
# Load filtered cluster abundance data
######################################

# load abundance data to repeat analysis
# The "f1" clusters have only had samples filtered out
# TARA prok srf will be analyzed here
load("~/Downloads/clstrs_tara_f2.RData")
#load("/bioinf/projects/megx/UNKNOWNS/Matt/clstrs_tara_f2.RData")
clstrs_tara_comp_stats <- clstrs_tara_comp_agg_stats_f2 %>%
  rename(label = sample_ID)

# Calculate the smallest sample library
minlib <- clstrs_tara_comp_stats %>%
  group_by(label) %>%
  summarise(N=sum(abun)) %>%
  .$N %>%
  min()

# Rounding function
myround <- function(x) { trunc(x + 0.5) }

# Round scale and round component abundances
clstrs_tara_comp_stats <- clstrs_tara_comp_stats %>%
  mutate(abun_scal = myround(proportion * minlib))

clstrs_tara_comp_stats$proportion %>% summary()


clstrs_eupc_stats <- clstrs_tara_comp_stats

# Filter for mean_prop > 2e-5
clstrs_eupc_stats_filtered <- clstrs_eupc_stats %>%
  filter(mean_proportion > 2e-5)

# Spread into matrix
clstrs_eupc_stats_filtered_df <- clstrs_eupc_stats_filtered %>%
  select(label, component, abun_scal) %>%
  #mutate(abun = as.integer(abun)) %>%
  spread(component, abun_scal, fill = 0) %>%
  as.data.frame()

# make label the rownames
rownames(clstrs_eupc_stats_filtered_df) <- clstrs_eupc_stats_filtered_df$label
clstrs_eupc_stats_filtered_df$label <- NULL

# Calculate Levins niche breadth (B) on samples
clstrs_eupc_final_filtered_B <- clstrs_eupc_stats_filtered %>%
  ungroup() %>%
  group_by(label) %>%
  mutate(proportion_n = (abun/sum(abun))) %>%
  ungroup() %>%
  group_by(component) %>%
  mutate(proportion = (abun/sum(abun)), n = length(unique(label))) %>%
  mutate(mean_proportion = mean(proportion_n)) %>%
  select(component, proportion, proportion_n, mean_proportion, n) %>%
  mutate(B = 1/(sum(proportion^2)), B_a = (B-1)/(n-1), B_a = ifelse(is.infinite(B_a), 0, B_a)) %>%
  #select(-prop ortion) %>%
  unique()

comm.tab <- clstrs_eupc_stats_filtered_df
comm.tab <- comm.tab[,which(colSums(comm.tab)>0)]
comm.tab <- comm.tab[which(rowSums(comm.tab)>0),]

####################################
# Spread data to make response table
####################################
# spread clstr_ID by abundance and fill "0" if no abundance
#clstrs_df <- dcast.data.table(sample_ID ~ component, fun.aggregate = "sum", fill = 0, data = as.data.table(clstrs), value.var = "abun_scal")

# retrieve rownames 
#rownames(clstrs_df) <- clstrs_df$sample_ID

#clstrs_df$sample_ID <- NULL

############################################################################
# Vegan null model package to calculate distribution of niche breadth "B" in 
# simulated population
############################################################################
#comm_tab <- clstrs_df

get_random <- function(x, y){
  null <- nullmodel(y, method = "quasiswap_count")
  res <- simulate(null, nsim=1)
  l <- as_data_frame(niche.width(res,  method = "levins"))
  return(l)
}

# levin.index.simul<-plyr::ldply(1:100, get_random, y = comm.tab, .parallel = T)
#levin.index.simul <- lapply(1, get_random, y = comm.tab) 
# Apply function 100x
levin.index.simul <- pbmclapply(1:100, get_random, y = comm.tab, mc.cores = 70) 

# bind the rows
levin.index.simul <- data.table::rbindlist(levin.index.simul)

# calculate the real niche breadth from the original response matrix
levin_index_real <- as.numeric(niche.width(comm_tab, method = "levins"))

# add colnames to from comm_tab to levin_index_simul
colnames(levin.index.simul) <- colnames(comm.tab)

levin.index.simul <- as.data.frame(levin.index.simul)

# Calculate mean
media <- apply(levin.index.simul, 2, mean)

# calculate the upper and lower quantiles
ci <- apply(levin.index.simul, 2, quantile, probs = c(0.025, 0.975))

# Make df
results <- data.frame(component = colnames(comm.tab), 
                      observed = levin.index.real, 
                      mean.simulated = media, 
                      lowCI = ci[1, ], 
                      uppCI = ci[2, ], 
                      sign = NA)

# resutls
results <- results %>% mutate(sign = case_when( observed > uppCI ~ 'Broad',
                                                observed < lowCI ~ 'Narrow',
                                                observed >= lowCI & observed <= uppCI ~ 'Non significant'))%>%
  tbl_df()




##############
# Save results
##############
save(results, file = "/bioinf/projects/megx/UNKNOWNS/Matt/TARA_srfprok_results_B.Rdata")

load("~/Downloads/niche_ALL_results_20180207.Rda")
#or
load("~/Dropbox/Public/matt_msc/data/niche_ALL_results_20180207.Rda")
results
dim(results)
str(results)

######
# Plot
######
# Separate into classes
results_k <- results %>%
  filter(grepl("k", component))
results_eu <- results %>%
  filter(grepl("eu", component))
results_gu <- results %>%
  filter(grepl("gu", component))

save(results, file =  "~/Dropbox/Public/matt_msc/data/nicheB_results_.Rdata")

# How many narrows and wides are in each category?
all_stats <- results %>% dplyr::select(component, sign) %>% unique %>% group_by(sign) %>% count() %>% mutate(category = "all")
k_stats <- results_k %>% dplyr::select(component, sign) %>% unique %>% group_by(sign) %>% count() %>% mutate(category = "k")
gu_stats <- results_gu %>% dplyr::select(component, sign) %>% unique %>% group_by(sign) %>% count() %>% mutate(category = "gu")
eu_stats <- results_eu %>% dplyr::select(component, sign) %>% unique %>% group_by(sign) %>% count() %>% mutate(category = "eu")

B_stats <- bind_rows(all_stats, k_stats, gu_stats, eu_stats)
B_stats$category <- factor(B_stats$category, levels = c("all", "k", "gu", "eu"))


results_k %>% dplyr::select(component, sign) %>% unique %>% .$sign %>% table() %>% prop.table()
results_gu %>% dplyr::select(component, sign) %>% unique %>% .$sign %>% table() %>% prop.table()
results_eu %>% dplyr::select(component, sign) %>% unique %>% .$sign %>% table() %>% prop.table()

# Plot proportions of Niche Breadth classification for each category
library(scales)
p0 <- ggplot(B_stats, aes(x = category, y = n, fill = sign)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom") +
  ggtitle("Niche breadth classification") +
  xlab("Component category") +
  ylab("Percent")

# add width and length!!!!!! 4:3 ratio length to width
ggsave(p0, filename = "../../../img/nicheB_proportions.png")

# plot all
all_plot <- ggplot(results %>% dplyr::select(mean_proportion, observed, sign) %>% unique, aes(mean_proportion, observed, fill = sign)) +
  geom_point(size=1.7, alpha = 0.9, shape=21, color = "#3A3A3A") + theme_light() +
  scale_x_log10() + 
  xlab("Mean proportion (log10)") + 
  ylab("Levin's niche breadth (B)") + 
  #ggtitle("Environmental unknowns PC distribution") + 
  guides(fill = guide_legend(override.aes = list(size=2), nrow = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position="top",
        legend.title=element_blank(), 
        legend.key=element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.4)),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.text.align=0) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6"))

eu_all_plot <- ggplot(results_eu %>% dplyr::select(mean_proportion, observed, sign) %>% unique, aes(mean_proportion, observed, fill = sign)) +
  geom_point(size=1.7, alpha = 0.9, shape=21, color = "#3A3A3A") + theme_light() +
  scale_x_log10() + 
  xlab("Mean proportion (log10)") + 
  ylab("Levin's niche breadth (B)") + 
  #ggtitle("Environmental unknowns PC distribution") + 
  guides(fill = guide_legend(override.aes = list(size=2), nrow = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position="top",
        legend.title=element_blank(), 
        legend.key=element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.4)),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.text.align=0) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6"))

gu_all_plot <- ggplot(results_gu %>% dplyr::select(mean_proportion, observed, sign) %>% unique, aes(mean_proportion, observed, fill = sign)) +
  geom_point(size=1.7, alpha = 0.9, shape=21, color = "#3A3A3A") + theme_light() +
  scale_x_log10() + 
  xlab("Mean proportion (log10)") + 
  ylab("Levin's niche breadth (B)") + 
  #ggtitle("Environmental unknowns PC distribution") + 
  guides(fill = guide_legend(override.aes = list(size=2), nrow = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position="top",
        legend.title=element_blank(), 
        legend.key=element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.4)),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.text.align=0) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6"))

k_all_plot <- ggplot(results_k %>% dplyr::select(mean_proportion, observed, sign) %>% unique, aes(mean_proportion, observed, fill = sign)) +
  geom_point(size=1.7, alpha = 0.9, shape=21, color = "#3A3A3A") + theme_light() +
  scale_x_log10() + 
  xlab("Mean proportion (log10)") + 
  ylab("Levin's niche breadth (B)") + 
  #ggtitle("Environmental unknowns PC distribution") + 
  guides(fill = guide_legend(override.aes = list(size=2), nrow = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position="top",
        legend.title=element_blank(), 
        legend.key=element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.4)),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.text.align=0) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6"))

p <- cowplot::plot_grid( all_plot + theme(legend.position="none") + ggtitle("All"),
                         k_all_plot + theme(legend.position="none") + ggtitle("Knowns"),
                         gu_all_plot + theme(legend.position="none") + ggtitle("Genomic unknowns"),
                         eu_all_plot + theme(legend.position="none") + ggtitle("Environmental unknowns"),
                         align = 'vh',
                         labels = c("A)", "B)", "C)", "D)"),
                         hjust = -1,
                         nrow = 1
)
legend_b <- cowplot::get_legend(all_plot + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p1 <- cowplot::plot_grid( p, legend_b, ncol = 1, rel_heights = c(1, .2))
p1


# add width and length!!!!!! 4:3 ratio length to width
ggsave(p1, filename = "../../../img/nicheB.png")

################################################
# Pick EUs and GUs from Narrow and Wide breadths
################################################
# gu
gu_narrow <- results_gu %>% 
  dplyr::select(-proportion, -proportion_n) %>%
  unique %>% 
  filter(mean_proportion > 1e-04, observed < 5)

clstrs_tara_comp_agg_stats_f2 %>%
  filter(component %in% gu_narrow$component) %>%
  group_by(sample_ID) %>%
  count() %>%
  ungroup() %>%
  top_n(n = 5, n) %>%
  arrange(desc(n))

gu_broad <- results_gu %>% 
  dplyr::select(-proportion, -proportion_n) %>%
  unique %>% 
  filter(mean_proportion > 2e-04, observed > 170)

clstrs_tara_comp_agg_stats_f2 %>%
  filter(component %in% gu_narrow$component) %>%
  group_by(sample_ID) %>%
  count() %>%
  ungroup() %>%
  top_n(n = 5, n) %>%
  arrange(desc(n))

eu_narrow <- results_eu %>% 
  dplyr::select(-proportion, -proportion_n) %>%
  unique %>% 
  filter(mean_proportion > 1e-04, observed < 5)

clstrs_tara_comp_agg_stats_f2 %>%
  filter(component %in% eu_narrow$component) %>%
  group_by(sample_ID) %>%
  count() %>%
  ungroup() %>%
  top_n(n = 5, n) %>%
  arrange(desc(n))
