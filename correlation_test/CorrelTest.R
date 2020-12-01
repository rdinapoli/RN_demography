library(rcarbon)
library(ggplot2)
library(gridExtra)
library(here)
load(here('R_imagefiles','variables.RData')) # Load Calibrated Dates and Custom Curves
source(here('src','randomThin.R'))

# Load data
x.norm=cal.combined.norm
x.nnorm=cal.combined.nnorm
bins = bins
ccurves = list(custom=customCurve)
timeRange = c(800,150)
palms = palm$PollenPerc
SOI = soi$SOIpr

# Create SPDs
set.seed(123)
index=randomThin(bins)
spd.norm =spd(x.norm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange, verbose=FALSE)
spd.nnorm =spd(x.nnorm[index],datenormalised=FALSE, spdnormalised = TRUE,timeRange=timeRange,verbose=FALSE)

# Create dataframe with spds and variables
nrm <- data.frame(x=spd.norm$grid$calBP, y=spd.norm$grid$PrDens, palms, SOI)
nnrm <- data.frame(x=spd.nnorm$grid$calBP, y=spd.nnorm$grid$PrDens, palms, SOI)

# Calculate correlation coefficients for entire time-series
spd.nnorm.palms.cor <- cor.test(nnrm$y, nnrm$palms, method="spearman")
spd.nnorm.palms.cor
spd.nnorm.soi.cor <- cor.test(nnrm$y, nnrm$SOI, method="spearman")
spd.nnorm.soi.cor
spd.norm.palms.cor <- cor.test(nrm$y, nrm$palms, method="spearman")
spd.norm.palms.cor
spd.norm.soi.cor <- cor.test(nrm$y, nrm$SOI, method="spearman")
spd.norm.soi.cor

# Extract estimates
nnorm.palms <- data.frame(name="nnorm.palms", SPD='non-normalized', r=spd.nnorm.palms.cor$estimate)#, low=spd.nnorm.palms.cor$conf.int[1], high=spd.nnorm.palms.cor$conf.int[2])
nnorm.soi <- data.frame(name="nnorm.soi", SPD='non-normalized', r=spd.nnorm.soi.cor$estimate)#, low=spd.nnorm.soi.cor$conf.int[1], high=spd.nnorm.soi.cor$conf.int[2])
norm.palms <- data.frame(name="norm.palms",SPD='normalized',r=spd.norm.palms.cor$estimate)#, low=spd.norm.palms.cor$conf.int[1], high=spd.norm.palms.cor$conf.int[2])
norm.soi <- data.frame(name='norm.soi',SPD='normalized',r=spd.norm.soi.cor$estimate)#, low=spd.norm.soi.cor$conf.int[1], high=spd.norm.soi.cor$conf.int[2])
cor.results <- rbind.data.frame(nnorm.palms, nnorm.soi, norm.palms, norm.soi)

# Calculation correlation coefficients for subsets
# Subset SPDs to pre- and post- 490 cal BP
nrm_pre <- subset(nrm, nrm$x >= 490)
nrm_post <- subset(nrm, nrm$x < 490)
nnrm_pre <- subset(nnrm, nnrm$x >= 490)
nnrm_post <- subset(nnrm, nnrm$x < 490)
# Calculate correlation for different time periods
spd.norm.palms.cor.pre <- cor.test(nrm_pre$x, nrm_pre$palms, method="spearman")
spd.norm.palms.cor.pre
spd.norm.palms.cor.post <- cor.test(nrm_post$y, nrm_post$palms, method="spearman")
spd.norm.palms.cor.post
spd.norm.soi.cor.pre <- cor.test(nrm_pre$y, nrm_pre$SOI, method="spearman")
spd.norm.soi.cor.pre
spd.norm.soi.cor.post <- cor.test(nrm_post$y, nrm_post$SOI, method="spearman")
spd.norm.soi.cor.post
spd.nnorm.palms.cor.pre <- cor.test(nnrm_pre$y, nnrm_pre$palms, method="spearman")
spd.nnorm.palms.cor.pre
spd.nnorm.palms.cor.post <- cor.test(nnrm_post$y, nnrm_post$palms, method="spearman")
spd.nnorm.palms.cor.post
spd.nnorm.soi.cor.pre <- cor.test(nnrm_pre$y, nnrm_pre$SOI, method="spearman")
spd.nnorm.soi.cor.pre
spd.nnorm.soi.cor.post <- cor.test(nnrm_post$y, nnrm_post$SOI, method="spearman")
spd.nnorm.soi.cor.post

# Create data.frames for plotting
norm.palms.pre <- data.frame(name='Pre-490 BP',SPD='normalized',r=spd.norm.palms.cor.pre$estimate)
norm.palms.post <- data.frame(name='Post-490 BP',SPD='normalized',r=spd.norm.palms.cor.post$estimate) 
norm.soi.pre <- data.frame(name='Pre-490 BP',SPD='normalized',r=spd.norm.soi.cor.pre$estimate)
norm.soi.post <- data.frame(name='Post-490 BP',SPD='normalized',r=spd.norm.soi.cor.post$estimate) 
nnorm.palms.pre <- data.frame(name='Pre-490 BP',SPD='non-normalized',r=spd.nnorm.palms.cor.pre$estimate) 
nnorm.palms.post <- data.frame(name='Post-490 BP',SPD='non-normalized',r=spd.nnorm.palms.cor.post$estimate) 
nnorm.soi.pre <- data.frame(name='Pre-490 BP',SPD='non-normalized',r=spd.nnorm.soi.cor.pre$estimate)
nnorm.soi.post <- data.frame(name='Post-490 BP',SPD='non-normalized',r=spd.nnorm.soi.cor.post$estimate) 
cor.results.prepost.palms <- rbind.data.frame(norm.palms.pre, norm.palms.post, nnorm.palms.pre, nnorm.palms.post)
cor.results.prepost.soi <- rbind.data.frame(norm.soi.pre, norm.soi.post, nnorm.soi.pre, nnorm.soi.post)

# plot spds by variables and correlation coefficients
p1 <- ggplot(data=nnrm, aes(x=x, y=y)) + 
  geom_line(aes(color=palms), size=2) + 
  scale_color_viridis_c() +
  labs(subtitle="a) non-normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=490, linetype='dashed', color='black', size=1)
p2 <- ggplot(data=nnrm, aes(x=x, y=y)) + 
  geom_line(aes(color=SOI), size=2) + 
  scale_color_viridis_c(option='magma')+
  labs(subtitle="b) non-normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=490, linetype='dashed', color='black', size=1)
#non-normalized SPD as a function of palms and SOI
p3 <- ggplot(data=nrm, aes(x=x, y=y)) + 
  geom_line(aes(color=palms), size=2) + 
  scale_color_viridis_c() +
  labs(subtitle="c) normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=490, linetype='dashed', color='black', size=1)
p4 <- ggplot(data=nrm, aes(x=x, y=y)) + 
  geom_line(aes(color=SOI), size=2) + 
  scale_color_viridis_c(option='magma')+
  labs(subtitle="d) normalized SPD") +
  xlab("cal BP")+
  ylab("SPD") +
  xlim(800,200) +
  geom_vline(xintercept=490, linetype='dashed', color='black', size=1)
cor.plot1 <- ggplot(data=cor.results.prepost.palms, aes(x=r, y=name, color=SPD)) + geom_point(size=2) +
  xlim(-1,1) +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  labs(x = "rho", y="Palm cover") +
  labs(subtitle="e) correlation between palm cover and SPDs")
cor.plot2 <- ggplot(data=cor.results.prepost.soi, aes(x=r, y=name, color=SPD)) + geom_point(size=2) +
  xlim(-1,1) +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  labs(x = "rho ", y="SOI index") +
  labs(subtitle="f) correlation between SOI and SPDs") 
grid.arrange(p1, p2, p3, p4, cor.plot1,cor.plot2, nrow=3)

jpeg(file=here('figures','temporary','correlation_results.jpeg'), width=12, height=8, units='in', res=300)
grid.arrange(p1, p2, p3, p4, cor.plot1,cor.plot2, nrow=3)
dev.off()

tiff(file=here('figures','temporary','correlation_results.tiff'), width=12, height=8, compression = 'lzw', units='in', res=600)
grid.arrange(p1, p2, p3, p4, cor.plot1,cor.plot2, nrow=3)
dev.off()
