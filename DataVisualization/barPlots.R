#df1 has 4 columns: Model, Pathway, Gain, Standard Deviation

ggplot(df1, aes(x=as.factor(Pathway), y=Gain, fill=Model)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=Gain-Std, ymax=Gain+Std), width=.2,position=position_dodge(.9)) +
  geom_hline(aes(yintercept=-5,linetype = "PSF2"),col='blue') +
  geom_hline(aes(yintercept=-15,linetype = "PSF5"),col='red')+
  scale_linetype_manual(name = "control", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c("blue", "red"))))+
  labs(title = "Change in Total Genetic Variance Across Models for Each Pathway", x = "Pathway",
       y = "Change in Variance")
