library(readxl)
library(tidyverse)

outlier <- function (df, value) {
  iqr <- boxplot(df[value], plot=FALSE)$stats
  df[which(df[value] >= iqr[5] | df[value] <= iqr[1]), value] = NA
  return(df[complete.cases(df),])
}

# 0 <= S <= 1 / S = sRGB value
linear_sRGB <- function (S) {
  if (S <= 0.04045){
    L <- S/12.92
  }
  else {
    L <- ((S+0.055)/1.055)^2.4
  }
  return (L)
}

path="C:\\Users\\user\\Desktop\\Paper\\Figure\\PRLib"
files <- paste0(path, "\\", list.files(path, "_index.xlsx"))

n=1
for (i in files) {
  assign(paste0("df", n), read_excel(i))
  n=n+1
}

tmp <- rbind(df1, df2)
for (i in 4:n-1) {
  tmp <- rbind(tmp, get(paste0("df", i)))
}

tmp[,c("Rmean", "Gmean", "R_R", "R_G")] <- linear_sRGB(tmp[,c("Rmean", "Gmean", "R_R", "R_G")]/255)
tmp["RFU"] <- ((tmp$Gmean+tmp$Rmean)*tmp$Area)/((tmp$R_R+tmp$R_G)*tmp$Rarea)

tmp["Full_Path"] <- apply(X = tmp["Full_Path"], 1,
                          FUN=function(i){
                            basename(i)
                            })

writexl::write_xlsx(tmp, paste0(path, "\\", "210826.xlsx"))

tmp %>% subset(second > 0.65 & count > 15 & score > 0.4) %>% ggplot() + 
  geom_point(aes(x=Rarea*R_R, y=RFU, color=Strain)) + 
  theme_bw() + facet_wrap(.~Part, ncol=5)

tmp %>% ggplot() + 
  geom_point(aes(x=Area, y=RFU, color=Strain)) + 
  theme_bw() + facet_wrap(.~Part, ncol=5)

tmp %>% subset(second > 0.65 & count > 15 & score > 0.4) %>% 
  group_by(Part, Strain) %>%
  summarise(n=n(), mean=mean(RFU), sd=sd(RFU)) %>%
  ggplot(aes(x=Part, y=mean)) + geom_bar(stat="identity") + facet_grid(.~Strain) +
  coord_flip() + geom_errorbar(
  aes(ymin=mean-sd, ymax=mean+sd),
  width=.3,
  size=.5,
  position=position_dodge(0.9))
