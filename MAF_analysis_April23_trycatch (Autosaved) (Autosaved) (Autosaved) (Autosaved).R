
setwd('~/Library/CloudStorage/OneDrive-MassGeneralBrigham/Crompton Folder For Active Files/Gabriela/ANBL1531/Pancancer_panel_May_2024/without_germline')
setwd('/Users/gv974/Library/CloudStorage/OneDrive-MassGeneralBrigham/Crompton Folder For Active Files/Gabriela/ANBL1531/Pancancer_panel_May_2024/with_germline')

list_files <- dir() #pattern='PZJ'
data = data.frame()
for (i in list_files) {
  tryCatch({
    # read the data
    file <- read.maf(i,vc_nonSyn=c('Frame_Shift_Ins','INS','DEL','Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','In_Frame_Ins','In_Frame_Del','Splice_Site',"3'UTR",'IGR',"5'UTR")) #,"3'UTR",'IGR',"5'UTR"
    file_data <- file@data
    maf_data <- read.maf(file_data,vc_nonSyn=c('Frame_Shift_Ins','INS','DEL','Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','In_Frame_Ins','In_Frame_Del','Splice_Site',"3'UTR",'IGR',"5'UTR")) #,"3'UTR",'IGR',"5'UTR"

    # create object name based on filename
    dataframe_name <- sub(".maf", "", i) #"\\-.*", "", i

    # name the object
    assign(dataframe_name, maf_data)
    maf_data <- maf_data@data[,c(1,5,6,7,9:13,16,37:42,80,81,82)] #141 c(1,5,6,7,9:13,16,37:42,80,81,82,141)
    data <- rbind(data,maf_data)
  }, error = function(e) {
    print(paste(i, "could not be read"))
  })
}

data$VAF <- data$t_alt_count/(data$t_alt_count + data$t_ref_count)*100
data$Tumor_Sample_Barcode<- gsub("^.*?_PDO-","",data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub("^.*?_P","P",data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub("\\_.*","",data$Tumor_Sample_Barcode)
data$time_point <- sub(".*-","",data$Tumor_Sample_Barcode)
data$sample <- sub("-.*","",data$Tumor_Sample_Barcode)
data$combined <- paste(data$Tumor_Sample_Barcode, data$Hugo_Symbol,data$Protein_Change, sep='_')

dim(data)
data.no_germline <- data
data.germline <- data

all_files <- dir()
exclude_pattern <- '(1)'
files_without_pattern <- all_files[!grepl(exclude_pattern, all_files)]


# mutation only in no_germline
data.difference <- data.no_germline[!data.no_germline$combined %in% data.germline$combined,]
data.difference2 <- data.difference[data.difference$TLOD < 50,]
dim(data.difference2)
# bind
data <- rbind(data.germline,data.difference2)
write_clip(data)

data$sample <- sub('-.*','',data$Tumor_Sample_Barcode)
data <- data[data$sample %in% pcp$V1,]

data$Time.point <- gsub(".*-","",data$Tumor_Sample_Barcode)su
-------


data$ichor <- ifelse(grepl("BL", data$Tumor_Sample_Barcode), 82,
                   ifelse(grepl("C2", data$Tumor_Sample_Barcode), 16,
                   ifelse(grepl("C3", data$Tumor_Sample_Barcode), 0,
                   ifelse(grepl("C4", data$Tumor_Sample_Barcode), 0,
                   ifelse(grepl("EOI", data$Tumor_Sample_Barcode), 0,
                   ifelse(grepl("SPC", data$Tumor_Sample_Barcode), 0,
                   ifelse(grepl("EPC", data$Tumor_Sample_Barcode), 0,
                   ifelse(grepl("RL", data$Tumor_Sample_Barcode), 72.37, NA))))))))

data$VAF <- data$t_alt_count/(data$t_alt_count + data$t_ref_count)*100
data$VAF_normalized <- 100/data$ichor*data$VAF


----

list_samples <- c('PAZCCH','PAZCFH','PAZELP','PAZFIW','PAZISC','PAZJFB','PAZJFJ','PAZMEX','PAZMNN','PAZNRG','PAZPJY','PAZPZJ','PAZRBA','PAZRTP','PAZSKT','PAZTBA','PAZTZE','PAZWVP')


data$Tumor_Sample_Barcode <- sub('\\-.*','',data$Tumor_Sample_Barcode)
data <- data[data$Tumor_Sample_Barcode %in% list_samples,]
summary(as.factor(data$Tumor_Sample_Barcode))
data$VAF <- data$t_alt_count/(data$t_alt_count + data$t_ref_count)*100

data <- data[data$VAF > 4,]





