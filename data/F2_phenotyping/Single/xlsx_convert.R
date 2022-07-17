path1="Z:\\09 Lab members\\박성군\\task\\Raw_data\\Image\\210630_DH5a_C2566\\DH5a_P50\\DH5a_P50_index.xlsx"
path2="Z:\\09 Lab members\\박성군\\task\\Raw_data\\Image\\210630_DH5a_C2566\\DH5a_P50\\DH5a_P50_index2.xlsx"

library(readxl)

df <- read_excel(path1, length(excel_sheets(path1)))
df2 <- read_excel(path2, length(excel_sheets(path2)))

colnames(df2)[1] <- "Index"

df[,c("Rmean", "Gmean", "R_R", "R_G")] = df2[df$Index,c("Rmean", "Gmean", "R_R", "R_G")]

new_index <- c("Index", "Full_Path", "X", "Y", "Area", "Rmean", "Gmean", "Rarea", "R_R", "R_G", "distance", "barcode", "tag_for", "tag_rev", "count", "Part", "score", "second", "Strain", "Date", "Time")
df <- df[,new_index]
writexl::write_xlsx(df, path1)