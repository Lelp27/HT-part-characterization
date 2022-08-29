# 220620 barplot

library(tidyverse)
library(readxl)
library(wesanderson)

data_path = "C://users//user//desktop//tcc_od.xlsx"

df <- read_excel(data_path, sheet=2)

df2 <- df %>% pivot_longer(cols=c('0h','8h', '16h'))


mycolor <- RColorBrewer::brewer.pal(5, 'Oranges')
mycolor[1] <- "#FFD1B7"

ggplot(df2, aes(x=name, y=value, fill=index, width=.6)) + geom_bar(position='stack', stat='identity') +
  scale_x_discrete(limits=c('16h', '8h', '0h')) + 
  coord_flip() + theme_bw() +
  theme(legend.position = 'none',
        axis.title=element_blank()) +
  scale_fill_manual(values = mycolor)
  #scale_fill_brewer(palette='Oranges')

df <- read_excel(data_path, sheet=3)
df <- df[c('Seq', 'Est')]

df %>% ggplot(aes(x=Seq, y=Est)) +
  geom_point(shape=24, fill= '#91f2ef', size=8) + theme_bw() +
  ylim(c(0,1)) + xlim(c(0,1)) + 
  geom_smooth(method='lm', se =FALSE, fullrange=TRUE, color='#2d2a2b') +
  theme(axis.title=element_blank(),
        axis.text=element_text(size=25))


cor(df['Seq'], df['Est'])
