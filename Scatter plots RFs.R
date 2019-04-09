#RFs for black bean by compound on a log scale
rfs <- read.csv("raw_data.csv", header=T)
bb <- read.csv("black_bean2.csv", header=T)
library(ggplot2)
library(grid)
library(scales)
source("Base breaks function.R")
theme_set(theme_bw(10))

bb <- rfs %>% filter(food == "Common bean [Black], whole, raw")

fig6a <- ggplot(bb, aes(x=mean.rf, y=reorder(compound, mean.rf), shape=process, colour=process)) +
  geom_point() + scale_shape_manual(values=c(5,6,17,19)) +
  geom_vline(xintercept=1, linetype="dashed") + theme_bw() +
  theme(legend.key.height=unit(4, "mm"), legend.title=element_blank(),
        legend.key=element_blank(), axis.title.y=element_blank(),
        #panel.border=element_rect(colour="black"), 
        panel.grid.major=element_blank()) +
  scale_x_continuous(name="Retention factor", trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  scale_colour_brewer(palette="Dark2") + #ggtitle("A")
  
ggsave("black bean slide.png", width=180, height=100, unit="mm", dpi=300, plot=fig6a)

#with lattice
library(latticeExtra)
plot <- xyplot(reorder(compound, mean.rf) ~ mean.rf, data=bb, groups=process,
               pch=c(1,2,3,4), xlab="Mean retention factor", ylab=NULL, auto.key=T,
               scales=list(x=list(log=10, equispaced.log=FALSE)))

trellis.device(device="png", filename="black bean slide lattice.png")
print(plot)
dev.off()

#--------------------------------------------------------------------------------------

#RFs of quercetin, Q3Rut, 5CQA and caffeic acid, on log scales
#Quercetin
q <- read.csv("q_only2.csv", header=T)
#q <- rfs %>% filter(compound == "Quercetin")
fig5a <- ggplot(q, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process, colour=process)) + 
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(axis.title=element_blank(), panel.border=element_rect(colour="black"), 
        legend.title=element_blank(), legend.key.height=unit(4, "mm"),
        panel.grid.major=element_blank(), legend.key=element_blank()) + #ggtitle("A") +
  scale_colour_brewer(palette="Dark2") +
  scale_shape_manual(values=c(21,16,15,17,14,6,5,19)) + scale_x_log10(breaks=c(0.2,0.5,1,2))

ggsave("Quercetin for slide.png", width=140, height=70, unit="mm", dpi=600, plot=fig5a)

#Q3Rut
q3r <- read.csv("q3r_only2.csv")
fig5b <- ggplot(q3r, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process, colour=process)) +
  geom_point() + 
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(4, "mm"),
        legend.title=element_blank(), #legend.position="none",
        panel.border=element_rect(colour="black"), legend.key=element_blank()) +
  scale_shape_manual(values=c(15,17,14,6,19)) + #ggtitle("B") +
  scale_x_log10(name="Retention factor", breaks=c(0.1,0.5,1,1.5)) +
  scale_colour_brewer(palette="Dark2") +
  geom_vline(xintercept=1, linetype="dashed")

ggsave("Q3Rut for slide.png", width=110, height=55, unit="mm", dpi=600, plot=fig5b)

#5-caffeoylquinic acid
fig4a <- ggplot(cqa, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process)) + 
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(legend.position="none", panel.border=element_rect(colour="black"),
        axis.title = element_blank(), panel.grid.major=element_blank(), legend.title=element_blank()) +
  scale_x_log10(breaks=c(0.05,0.1,0.2,0.5,1,2)) + ggtitle("A") +
  scale_shape_manual(values=c(21,15,17,14,18,6,5,13,19))

ggsave("fig4a revised.png", width=80, height=60, unit="mm", dpi=600, plot=fig4a)

#caffeic acid RFs on a log scale
ca <- read.csv("ca_only2.csv", header=T)
fig4b <- ggplot(ca, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process)) +
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(legend.key.height=unit(4, "mm"), panel.border=element_rect(colour="black"), 
        axis.title.y = element_blank(), legend.title=element_blank(),
        panel.grid.major=element_blank(), legend.key=element_blank()) +
  scale_x_log10(name="Retention factor", breaks=c(0.2,1,0.5,1,2,5)) +
  scale_shape_manual(values=c(21,16,15,17,14,18,6,5,13,19)) + ggtitle("B")

ggsave("fig4b revised.png", width=120, height=70, unit="mm", dpi=600, plot=fig4b)

#---------------------------------------------------------------------------------------

#potato scatter
potato <- read.csv("potato.csv", header=T)
#may need to find location of file to source
source("base_breaks.R")

fig6b <- ggplot(potato, aes(x=mean.rf, y=reorder(compound, mean.rf), shape=process)) +
  geom_point() + scale_shape_manual(values=c(21,5,15,18,19)) + 
  geom_vline(xintercept=1, linetype="dashed") +
  ggtitle("B") +
  theme(legend.key.height=unit(4, "mm"), axis.title.y=element_blank(), 
        panel.border=element_rect(colour="black"), legend.title=element_blank(),
        panel.grid.major=element_blank(), legend.key=element_blank()) +
  scale_x_continuous(name="Retention factor", trans=log_trans(),
                     breaks=base_breaks(), labels=prettyNum)

ggsave("fig6b revised.png", width=171, height=50, unit="mm", dpi=600, plot=fig6b)

#--------------------------------------------------------------------------------------------

#facet scatter RFs folin. subset Folin data only
rfs.folin <- subset(rfs, method == "Folin assay")
rfs.folin1 <- rfs %>% filter(method == "Folin assay", process == "Boiled")
ggplot(rfs.folin, aes(x=process, y=mean.rf)) + geom_point() + facet_grid(food.group ~ .) + theme_bw() +
  geom_jitter(position=position_jitter(w=0.2, h=0.15)) + scale_y_continuous(name="Mean retention factor") +
  scale_x_discrete(limits=c("Baked", "Blanched", "Boiled", "Fried", "Grilled",
                            "Microw.", "Pr. boiled", "Pr. steamed", "Steamed"))

#with lattice
library(latticeExtra)
xyplot(process ~ mean.rf, groups=food.group, data=rfs.folin)

#--------------------------------------------------------------------------------------------

#scatter RFs by process
ggplot(rfs, aes(x=process, y=mean.rf, colour=compound.class)) +
  geom_jitter(position=position_jitter(w=0.2, h=0.15)) + 
  #geom_hline(yintercept=1, linetype="dashed") + 
  #annotate("text", x=0.95, y=1.2, label="Mean RF=1") +
  theme_bw() + labs(x="Process", y="Aggregated retention factor") + 
  scale_x_discrete(limits=c("Baked", "Blanched", "Boiled", "Fried", "Grilled",
                            "Microw.", "Pr. boiled", "Pr. steamed", "Steamed")) +
  scale_y_log10(name="Mean retention factor") + theme_bw() +
  theme(legend.position=c(0.13,0.8))

ggsave("scatter rfs by process.png", width=6.5, height=4, unit="in", plot=processes)