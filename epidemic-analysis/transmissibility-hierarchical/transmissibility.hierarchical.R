library(stringr)
library(readr)
library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(dplyr)
library(tidyr)
library(tidytree)
library(openxlsx)
library(data.table)
library(magrittr)
library(tidyverse)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)

###inputs###
stan_f.name <- file.path(cmdstan_path(), 'data/multinomial_time_series_hierarchical.formatted.stan') #should-be-done-separately

#model
multi_nomial_model <- cmdstan_model(stan_f.name)

#general
core.num <- 4
pango.reference <- "EG.5.1"

#period to be analyzed
date.end <- as.Date("2023-10-04")
date.start <- as.Date("2023-04-01")

#Transmissibility
bin.size <- 1
generation_time <- 2.1

final.metadata.filtered <- fread('final.metadata.filtered.tsv', header=T, sep="\t", quote="",check.names=T)
metadata.BA.2.86 <- final.metadata.filtered %>% filter(clade_nextstrain %in% c("BA.2.86"))
count.country.BA.2.86.df <- metadata.BA.2.86 %>% group_by(country) %>% summarize(count = n()) %>% arrange(desc(count))
count.country.BA.2.86.df.filtered <- count.country.BA.2.86.df %>% filter(count >= 20) %>% arrange(desc(count))
count.country.BA.2.86.df.filtered <- count.country.BA.2.86.df.filtered[-3,] #removing-South-Africa
filtered.countries.BA.2.86 <- count.country.BA.2.86.df.filtered$country
filtered.countries.all <- final.metadata.filtered %>% filter(country %in% filtered.countries.BA.2.86)
strain.counts <- filtered.countries.all %>% group_by(country, group) %>% summarize(count = n()) %>% arrange(country, desc(count))
subset.strain.counts <- strain.counts %>% pivot_wider(names_from = group, values_from = count, values_fill = 0)

pango.interest.v <- c("XBB.1.5","XBB.1.16","EG.5.1","BA.2.86") #remove XBB.1
filtered.countries.interest <- filtered.countries.all %>% filter(group %in% pango.interest.v)
filtered.countries.interest <- filtered.countries.interest %>% mutate(group = factor(group,levels=pango.interest.v))
filtered.countries.interest <- filtered.countries.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
filtered.countries.interest <- filtered.countries.interest %>% filter(!is.na(date.bin))

out.name <- 'metadata.used_for_analysis.hierarchical.txt'
write.table(filtered.countries.interest,out.name,col.names=T,row.names=F,sep="\t",quote=F)

metadata.filtered.interest.bin <- filtered.countries.interest %>% group_by(date.bin.num,group,country) %>% summarize(count = n()) %>% ungroup()
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=group,value = count)
metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread %>% mutate(country = factor(country))
country <- as.numeric(metadata.filtered.interest.bin.spread$country)

X <- metadata.filtered.interest.bin.spread$date.bin.num

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num,-country)
Y <- Y[,c(pango.reference,colnames(Y)[-which(colnames(Y)==pango.reference)])]

count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()

group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))
country.df <- data.frame(country_Id = 1:length(levels(metadata.filtered.interest.bin.spread$country)),
                         country = levels(metadata.filtered.interest.bin.spread$country))


Y_sum.v <- apply(Y,1,sum)


###fitting###
data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  projection = nrow(Y),
                  D = max(country),
                  X = X,
                  Y = Y,
                  country = country,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = c(Y_sum.v))

fit.stan <- multi_nomial_model$sample(
  data=data.stan,
  iter_sampling=2000,
  iter_warmup=1000,
  seed=1234,
  parallel_chains = 4,
  adapt_delta = 0.99,
  max_treedepth = 15,
  #pars=c('b_raw'),
  chains=4)

##########outputs##########
#growth rate mean
stat.info.mean <- fit.stan$summary("growth_rate_mean") %>% as.data.frame()
stat.info.mean.q <- fit.stan$summary("growth_rate_mean", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame()
stat.info.mean <- merge(stat.info.mean,stat.info.mean.q,by="variable")
stat.info.mean <- stat.info.mean %>% mutate(group_Id = str_match(variable,'growth_rate_mean\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)

stat.info.mean.merged <- merge(stat.info.mean,group.df,by="group_Id") %>% select(-group_Id,-variable)
stat.info.mean.merged <- stat.info.mean.merged %>% arrange(desc(mean))


#growth_rate each country
stat.info.each <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info.each.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame()
stat.info.each <- merge(stat.info.each,stat.info.each.q,by="variable")

stat.info.each <- stat.info.each %>% mutate(country_Id = str_match(variable,'growth_rate\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(), group_Id = str_match(variable,'growth_rate\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric() + 1)
stat.info.each.merged <- stat.info.each %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id") %>% select(-group_Id,-country_Id,-variable)


out.name <- paste('data/Re-data/growth_rate.mean.hierarchical.txt',sep="")
write.table(stat.info.mean.merged,out.name,col.names=T,row.names=F,sep="\t",quote=F)

out.name <- paste('data/Re-data/growth_rate.each_country.hierarchical.txt',sep="")
write.table(stat.info.each.merged,out.name,col.names=T,row.names=F,sep="\t",quote=F)

#color-matching
color.matching <- c(
  "outgroup" = "grey",
  "Alpha" = "#1C67A8",
  "Beta" = "#39A04A",
  "Gamma" = "#D7C046",
  "Delta" = "#E77A25",
  "BA.1" = "#339680",
  "BA.2" = "#02686B",
  "BA.5" = "#803345",
  "BA.2.75" = "#898731",
  "BQ.1" = "#67468C",
  "XBB.1" = "#736357",
  "XBB.1.5" = "#981D2A",
  "XBB.1.16" = "#D95A24",
  "EG.5.1" = "#DE0303",
  "BA.2.86" = "#1D9FF0",
  "BA.4" = "#8C2788"
)

###growth rate_mean
#growth rate
draw.df.growth_rate <- fit.stan$draws("growth_rate_mean", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate_mean\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% group_by(group) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))
draw.df.growth_rate.long <- rbind(data.frame(group=pango.reference,value=1),draw.df.growth_rate.long)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group = factor(group,levels=pango.interest.v))

col.v <- c("XBB.1.5" = "#981D2A", "XBB.1.16" = "#D95A24", "EG.5.1" = "#DE0303","BA.2.86" = "#1D9FF0")
draw.df.growth_rate.long <- draw.df.growth_rate.long
g <- ggplot(draw.df.growth_rate.long,aes(x=group,y=value,color=group,fill=group))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.6,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1.5)
g <- g + scale_color_manual(values=col.v)
g <- g + scale_fill_manual(values=col.v)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + xlab('') + ylab(paste('Relative Re (',pango.reference,')',sep=""))
g <- g + theme(legend.position = 'none')
g <- g + scale_y_continuous(limits=c(0.8, 1.2), breaks=seq(0.8, 1.2, by=0.2))
g

pdf.name <- 'data/Re-data/growth_rate.global_mean.pdf'

pdf(pdf.name,width=1.8,height=3.5)
plot(g)
dev.off()

###growth rate for each country
draw.df.growth_rate.each <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.each.long <- draw.df.growth_rate.each %>% gather(key = class, value = value)

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% mutate(country_Id = str_match(class,'growth_rate\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(), group_Id = str_match(class,'growth_rate\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id") %>% select(-group_Id,-country_Id,-class)
draw.df.growth_rate.each.long <- rbind(data.frame(value=1,group=pango.reference,country=country.df$country),draw.df.growth_rate.each.long)

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long
draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% group_by(group,country) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% mutate(group = factor(group,levels=pango.interest.v))


g <- ggplot(draw.df.growth_rate.each.long,aes(x=group,y=value,color=group,fill=group))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.4,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1)
g <- g + scale_color_manual(values=col.v)
g <- g + scale_fill_manual(values=col.v)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + xlab('') + ylab('Relative Re (EG.5.1)')
g <- g + theme(legend.position = 'none')
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + scale_y_continuous(limits=c(0.8, 1.2), breaks=seq(0.8, 1.2, by=0.2))

g <- g + facet_wrap(~country,ncol=7)
g

pdf.name <- 'data/Re-data/reproduction_number.each_country.pdf'

pdf(pdf.name,width=6,height=5)
plot(g)
dev.off()


#theta
data_Id.df <- data.frame(data_Id = 1:length(X), country_Id = country, date_Id = X, Y_sum = Y_sum.v, date = as.Date(X,origin=date.start))


data.freq <- metadata.filtered.interest.bin %>% rename(group = group) %>% filter(group != "others") %>% group_by(country,date.bin.num) %>% mutate(freq = count / sum(count))
data.freq <- merge(data.freq,data_Id.df,by.x="date.bin.num",by.y="date_Id")


draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.theta.long <- draw.df.theta %>% gather(key = class, value = value)
draw.df.theta.long <- draw.df.theta.long %>% mutate(data_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())


draw.df.theta.long <- draw.df.theta.long %>% inner_join(data_Id.df,by="data_Id")
draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(group_Id, country_Id, date) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))

draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id")
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(data.freq %>% select(group,count,freq,date,country),by=c("date","group","country"))

draw.df.theta.long.sum.filtered <- draw.df.theta.long.sum %>% mutate(group = factor(group,levels=pango.interest.v))

g <- ggplot(draw.df.theta.long.sum.filtered,aes(x=date, y = mean, fill=group, color = group))
g <- g + geom_ribbon(aes(ymin=ymin,ymax=ymax), color=NA,alpha=0.4)
g <- g + geom_line(size=0.3)
g <- g + scale_x_date(date_labels = "%m", date_breaks = "1 months", date_minor_breaks = "1 month")
g <- g + theme_set(theme_classic(base_size = 6, base_family = "Helvetica")) #over-all
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8) #country-names
)
g <- g + facet_wrap(~country,ncol=7)
g <- g + scale_color_manual(values = col.v)
g <- g + scale_fill_manual(values = col.v)
g <- g + scale_size_continuous(range = c(0.2, 4))
g

pdf.name <- 'data/Re-data/theta.each_country.pdf'
pdf(pdf.name,width=10,height=2)
plot(g)
dev.off()

#done
