#Set the working directory
setwd("/Volumes/NGS-data/Cole/20200213")
#Import necessary libraries
library("dplyr")
library("plyr")
library(Biostrings)
library(ggplot2)
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(ShortRead)
library(ggseqlogo)

####################################################################################################
####################################################################################################
#find_CDR, calculate lenght of each CDR， and plot the sequence logo based on differnt length########
####################################################################################################
####################################################################################################

setwd("/Volumes/NGS-data/Cole/20200213")
#Read file names
raw <- list.files(pattern="\\.fastq.gz$")

count.summary <- matrix(ncol= 3) %>% as.data.frame()
#Primary for loop
for (n in 1:length(raw)) {
  #Read in Fastq file
  cat(paste("Currently reading in",raw[n],"\n"))
  write.table(t(c("File name",raw[n])), file = paste("processing_files/",strtrim(raw[n],13),"_","Processing.csv"), append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
  fq <- readFastq(raw[n])
  
  #Convert fastq file to sequences
  sequences <- sread(fq)
  cat(paste("There are",length(sequences),"reads in",raw[n],"\n"))
  write.table(t(c("Total counts",length(sequences))), file = paste("processing_files/",strtrim(raw[n],13),"_","Processing.csv"), append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
  
  #Clear the memory
  rm(fq)
  gc()
  
  #Tabulate the data
  tbl <- as.data.frame(table(sequences))
  cat(paste("Tabulation complete. There are", nrow(tbl),"unique sequences. That's", (nrow(tbl)/length(sequences))*100,"percent of the total reads.","\n"))
  write.table(t(c("Unique counts",nrow(tbl))), file = paste("processing_files/",strtrim(raw[n],13),"_","Processing.csv"), append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
  rm(sequences)
  gc()
  
  #Order by frequency and write to csv
  ord_tbl <- tbl[order(-tbl[,2]),]
  write.csv(ord_tbl,paste("ordered_tables/",strtrim(raw[n],13),".csv",sep=""))
  
  #Finishing time
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste("Processing",raw[n],"took",as.numeric(time.taken),"seconds.","\n"))
  write.table(t(c("Processing time",time.taken)), file = paste("processing_files/",strtrim(raw[n],11),"_","Processing.csv"), append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
  write.table("", file = paste("processing_files/",strtrim(raw[n],13),"_","Processing.csv"), row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
  
  #Clean-up (except raw)
  rm(tbl,ord_tbl)
  gc()
  
  #End the loop
}

######################
##### REFERENCES #####
######################

#https://stackoverflow.com/questions/16566799/change-variable-name-in-for-loop-using-r
#https://www.rdocumentation.org/packages/ShortRead/versions/1.30.0/topics/FastqFile-class
#https://stackoverflow.com/questions/9107695/r-comment-out-block-of-code





#list names of all files
file_name <- list.files(pattern = "\\.csv$")

read <- matrix(nrow = length(file_name), ncol = 5) %>% as.data.frame()
colnames(read) <- c("ID", "Total", "No_error", "No_stop", "In_frame")

for (i in 1:length(file_name)){
  
  NT_seq<-read.csv(file_name[i], stringsAsFactors = F)
  #ID
  read$ID[i] <- file_name[i]
  #total read
  read$Total[i] <- nrow(NT_seq)
  #remove reading error
  NT_seq <- NT_seq[!grepl(NT_seq$sequences, pattern="N"),]
  #No_error
  read$No_error[i] <- nrow(NT_seq)
  #translate
  NT_seq$AA_seq<- translate(DNAStringSet(NT_seq$sequences))
  #remove any stop_codon* 
  NT_seq <-subset(NT_seq, grepl("[*]", NT_seq$AA_seq)==FALSE)
  read$No_stop[i] <- nrow(NT_seq)
  # in_frame
  NT_seq <- NT_seq[grepl(NT_seq$AA_seq, pattern="DYWGQGT"), ]
  read$In_frame[i] <- nrow(NT_seq)
  #get the L3 or H3, please change the flanking sequences
  NT_seq$H3<- sub(as.character(NT_seq$AA_seq), pattern="DYWGQGT\\S*", replacement = "")
  NT_seq$length <- nchar(NT_seq$H3)
  NT_seq <- NT_seq[, c("sequences","AA_seq","Freq", "H3","length")]
  
  #NT_seq <- NT_seq[NT_seq$Freq>=5,]
  write.csv(NT_seq, file_name[i])
  
}

write.csv(read, "NGS_counts_summary.csv")




##########################################################
##########################################################
#alignment according to length and calcualt the length distribution
##########################################################
##########################################################


len_sum <- as.data.frame(matrix(nrow = 1, ncol = 18))
colnames(len_sum) <- 1:18

for (i in 1:length(file_name)){
  a<-read.csv(file_name[i], stringsAsFactors = F)
  len <- unlist(unique(a$length))
  for (j in 1:length(len)){
    b <- a[a$length==len[j],]
    len_sum[i, len[j]] <- nrow(b)
    
    ggseqlogo(b$H3, method= "bits")+
      labs(title=paste0(file_name[i], "_", len[j]), x="\nPosition", y = "Bits\n") +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.title  = element_text(size=10, face="bold", hjust = 0.1),
            legend.text = element_text(size=9, face="bold"), 
            plot.title= element_text(hjust = 0.5, face ="bold", size = 12), 
            axis.text=element_text(size=11,face="bold"),
            axis.title=element_text(size=12,face="bold")
      )
    ggsave(paste0(file_name[i], "_", len[j],".tiff"), plot = last_plot(), device = "tiff", path = NULL,
           scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
           dpi = 600, limitsize = TRUE)
  }
  
}


rownames(len_sum) <- file_name
len_sum <- len_sum[5:19]
write.csv(len_sum, "CDR_length_summary.csv")





#############################
#############################
#bar plot of length distribution
#############################
#############################

c <- read.csv("CDR_length_summary.csv")
colnames(c) <- sub(colnames(c), pattern = "X",replacement = "")
c <- c[, c(1, 5:19)]
c.m <- melt(c)
#colnames(c.m) <- c("ID", "variable", "value")


c.m$value <- c.m$value/1000000
ggplot(c.m, aes(x=variable, y=value)) +
  geom_bar(stat="identity", color="black", fill="steelblue", position=position_dodge(), width = .7)+
  theme_minimal()+
  labs(title="", x="\nH3 Length", y = "Counts\n") +
  theme(
        legend.title  = element_text(size=10, face="bold", hjust = 0.1),
        legend.text = element_text(size=9, face="bold"), 
        plot.title= element_text(hjust = 0.5, face ="bold", size = 12), 
        axis.text=element_text(size=11,face="bold"),
        axis.title=element_text(size=12,face="bold"))#+
  scale_y_continuous(breaks=seq(-20,20,4), limits = c(0, ))

ggsave("L3_lenth_distribution.tiff", plot = last_plot(), device = "tiff", path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)  





##############################
#calculate theoretical diversity
##############################

x <- read.csv("/Volumes/NGS-data/Cole/20200213/CDR_length_summary.csv")
for (i in 3:ncol(x)){
  x[2,i] <- 11^(i-2)*8
}

write.csv(x, "H3_length_summary.csv")

x <- read.csv("H3_length_summary.csv")
x <- read.csv("/Volumes/NGS-data/Cole/20200213/CDR_length_summary.csv", stringsAsFactors = F)
x[2,1] <- NA
x <- t(x) %>% as.data.frame()
colnames(x) <- c("Counts", "position")
x <- x[5:19,]
x$position <- seq(4,18,1)

x$NGS <- sprintf(x$NGS, fmt = '%#.1f')
#x$position <- paste("AA", x$position)

x$counts <- as.numeric(unlist(x$Counts))


ggplot(c.m, aes(x=variable, y=value)) +
  geom_bar(stat="identity", color="white", fill="steelblue", position=position_dodge(), width = .7)+
  theme_minimal()+
  labs(title="H3 length distribution", x="\nH3 Length", y = "Counts (x 10E6)\n") +
  theme(
    legend.title  = element_text(size=25, face="bold", hjust = 0.1),
    legend.text = element_text(size=20, face="bold"), 
    plot.title= element_text(hjust = 0.5, face ="bold", size = 12), 
    axis.text=element_text(size=18,face="bold"),
    axis.title=element_text(size=20,face="bold"))#+
   #scale_x_discrete(limits=c(4:18))


ggsave("L3_lenth_distribution_10E6.tiff", plot = last_plot(), device = "tiff", path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 600, limitsize = TRUE)  
#############################
#############################
###########heatmap###########
#############################
#############################

setwd("/Users/jiezhou/Box/Heatmap")
file_name <- list.files(pattern = ".csv")

for (k in 2:length(file_name)){
  trial <- read.csv(file_name[k])
  #trial$RPM <- trial$Freq/(sum(trial$Freq)/1000000)
  #add a column of RPM
  # write.csv(trial, file_name[k])
  #trial <- trial[trial$RPM>2, ]
  #column name, and define the sequence
  
  AA<-c("A", "C", "D", "E", "F", "G", "H","I",  "K", "L", "M", "N", "P", "Q", "R", "S", "T","V", "W", "Y")
  # a summarized table
  name <- sub(file_name[k], pattern = "H3_", replacement = "")
  name <- sub(name, pattern = ".csv", replacement = "") %>% as.numeric()
  
  
  calculation<- matrix(nrow=name, ncol=20) %>%as.data.frame
  colnames(calculation)<-AA
  # summarized how many each AA in each position
  for (i in 1:name) {
    temp<-substr(trial$H3, i, i)
    uniqueaa<- temp %>%  unique() %>%data.frame(stringsAsFactors = F)
    colnames(uniqueaa) <- "uniqueaa"
    result<-apply(uniqueaa,1,function(x) str_count(temp, x["uniqueaa"]))
    colnames(result) <- uniqueaa$uniqueaa
    extra <- AA[!AA%in%uniqueaa$uniqueaa]
    if (!length(extra)==0){
      extra1 <- result[,1:(length(extra)+1)]
      extra1[extra1!=0] <- 0
      colnames(extra1 ) <- c(extra, "temp")
      result <- cbind(result,extra1)
      result <- result[, 1:20]
    }
    result<-result[, AA]
    for (j in 1:20) {
      calculation[i,j]<-sum(result[,j])
    }
  }
  calculation$total<-sum(calculation[, 1:20])/name
  calculation_percentage<- calculation/calculation$total[1]*100
  
  calculation_percentage <- calculation_percentage[, 1:20] %>% t() %>% as.data.frame()
  Zscore <- calculation_percentage
  
  for (x in 1:name){
    for (y in 1: 20) {
      Zscore[y,x]<- (calculation_percentage[y,x]-5)/sd(calculation_percentage[, x])
    }
  }
  
  
  Zscore$AA <- rownames(Zscore)
  Zscore.m <- melt(Zscore)
  
  Zscore.m$AA <- factor(Zscore.m$AA, levels = rev(c("A", "C", "D", "E", "F", "G", "H","I",  "K", "L", "M", "N", "P", "Q", "R", "S", "T","V", "W", "Y")))
  ggplot(Zscore.m, aes(variable, AA, fill=value)) + 
    geom_tile( color = "white" , size=1)+
    scale_fill_gradient2( low = "blue4", mid = "white",
                          high = "lightcoral", space = "Lab",
                          na.value = "grey50", name="Z-Score")+
    #theme_minimal()+
    labs(title="", x="", y = "")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
          legend.title  = element_text(size=20, face="bold", hjust = 0.1),
          legend.text = element_text(size=13, face="bold"), 
          plot.title= element_text(hjust = 0.5, face ="bold", size = 20), 
          axis.text=element_text(size=18, face="bold"),
          axis.text.x = element_text(angle = 90))+
    #guides(fill=guide_legend(title="Z-score"))+
    coord_fixed(ratio = 1)
  
  ggsave(paste0("heatmap_",name, ".tiff"), plot = last_plot(), device = "tiff", path = NULL,
         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
         dpi = 300, limitsize = TRUE)
  
}





#################################
##percentage of each AA
#################################

trial <- read.csv("/Users/jiezhou/H3_18.csv", stringsAsFactors = F)

AA<-c("A", "C", "D", "E", "F", "G", "H","I",  "K", "L", "M", "N", "P", "Q", "R", "S", "T","V", "W", "Y")

calculation<- matrix(nrow=18, ncol=20) %>%as.data.frame
colnames(calculation)<-AA

for (i in 1:18) {
  temp<-substr(trial$H3, i, i)
  uniqueaa<- temp %>%  unique() %>%data.frame(stringsAsFactors = F)
  colnames(uniqueaa) <- "uniqueaa"
  result<-apply(uniqueaa,1,function(x) str_count(temp, x["uniqueaa"]))
  colnames(result) <- uniqueaa$uniqueaa
  extra <- AA[!AA%in%uniqueaa$uniqueaa]
  if (!length(extra)==0){
    extra1 <- result[,1:(length(extra)+1)]
    extra1[extra1!=0] <- 0
    colnames(extra1 ) <- c(extra, "temp")
    result <- cbind(result,extra1)
    result <- result[, 1:20]
  }
  result<-result[, AA]
  for (j in 1:20) {
    calculation[i,j]<-sum(result[,j])
  }
}
calculation$total<-sum(calculation[, 1:20])/18
calculation_percentage<- calculation/calculation$total[1]*100

calculation_percentage <- calculation_percentage[, 1:20] %>% t() %>% as.data.frame()
for (i in 1: 18){
  calculation_percentage$AA_perc[i] <- sum(calculation_percentage[i,1:18])
}

write.csv(calculation,"aa.csv")



p <- read.csv("/Volumes/NGS-data/Cole/20200213/AA_perce.csv", stringsAsFactors = F)
p <- p[,1:4]

ggplot(p,aes(x=Experiment, y=Theoretical) )+
  geom_point(aes(colour = Position), size=4)+
  # geom_text(label = x$X, nudge_y=-1)+
  #annotate("text", label=lm_eqn(x[,3:4]), parse=TRUE, x=Inf, y=Inf, hjust=0.5, vjust=1)+
  #theme_minimal()+
  #theme_classic()+
  labs(title="", y="Theoretical (%) \n", x = "\nExperimental (%)")+
  theme(legend.title  = element_text(size=18, face="bold", hjust = 0.5),
        legend.text = element_text(size=15, face="bold"),     
        plot.title= element_text(hjust = 0.5, face ="bold", size = 15),
        axis.text=element_text(size=23, face="bold"),         
        axis.title=element_text(size=27,face="bold")
  )+
  geom_text(label = p$AA, nudge_y=-1.5, nudge_x=1.5,size =8, face="bold")+
  geom_smooth(method = lm, se=T, color="black", size=1)+
  scale_y_continuous(limit=c(0,55),breaks=seq(0,100,10))+
  scale_x_continuous(limit=c(0,55),breaks=seq(0,100,10))+
  coord_fixed(ratio = 1)


ggsave("exp_ther_plot6.tiff", plot = last_plot(), device = "tiff", path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)














