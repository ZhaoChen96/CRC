killDbConnections = function () {
  all_cons = dbListConnections(MySQL())
  for(con in all_cons){dbDisconnect(con)}
}
args=commandArgs(T)
signature = args[1]       ## signature genelist
dataset = args[2]         ## selected datasets
parameter = args[3]       ## parameters
signatures = strsplit(signature,",")[[1]]       ## split the signature genelist into array
if(is.na(signatures[4])){signatures[4] = ""}
datasets = strsplit(dataset,",")[[1]]           ## split the datasets info into array
datasets = datasets[datasets != ""]
parameters = strsplit(parameter,",")[[1]]
dbs = "GE_SF"             ## database
symbol1 = parameters[1]
symbol2 = parameters[2]
symbol3 = parameters[3]
symbol4 = parameters[4]
method = parameters[5]          ## "pearson" or "spearman" or "kendall"
outputdir = parameters[6]       ## outputdir


.libPaths(c("--------",.libPaths()))      ## add the specific directory of RMySQL R package
suppressPackageStartupMessages(library("RMySQL"))
killDbConnections()
#mydb = dbConnect(MySQL(), user='--------', password='--------', dbname=dbs)     ## connet to mysql database
mydb = dbConnect(MySQL(), user='--------', password='--------', dbname=dbs,unix.socket="--------") 

## Gene A:
if(signatures[3] != ""){
  symbol1 = paste(symbol1,symbol3,sep = "/")
  signatures_a = c(signatures[1],signatures[3])
}else{signatures_a = signatures[1]}
df_a = as.matrix(array(data = 0,dim = c(0,length(signatures_a))))
for(i in 1:length(datasets)){
  table = datasets[i]
  df_t=t(dbGetQuery(mydb,paste("SELECT * FROM ",table," WHERE geneid IN ('",paste(signatures_a, collapse = "','"),"')",sep="")))
  colnames(df_t) = df_t[1,]
  df_t = df_t[-1,,drop = F]
  df_a = rbind(df_a,df_t)
}
storage.mode(df_a) = "numeric"
df_a = df_a[,signatures_a,drop = F]
if(signatures[3] != ""){
  df_a = log2(df_a + 0.001)
  df_a = df_a[,1,drop = F] - df_a[,2,drop = F]
  df_a = 2^df_a
  colnames(df_a) = signatures_a[1]
}
## Gene B:
if(signatures[4] != ""){
  symbol2 = paste(symbol2,symbol4,sep = "/")
  signatures_b = c(signatures[2],signatures[4])
}else{signatures_b = signatures[2]}
df_b = as.matrix(array(data = 0,dim = c(0,length(signatures_b))))
for(i in 1:length(datasets)){
  table = datasets[i]
  df_t=t(dbGetQuery(mydb,paste("SELECT * FROM ",table," WHERE geneid IN ('",paste(signatures_b, collapse = "','"),"')",sep="")))
  colnames(df_t) = df_t[1,]
  df_t = df_t[-1,,drop = F]
  df_b = rbind(df_b,df_t)
}
storage.mode(df_b) = "numeric"
df_b = df_b[,signatures_b,drop = F]
if(signatures[4] != "" ){
  df_b = log2(df_b + 0.001)
  df_b = df_b[,1,drop = F] - df_b[,2,drop = F]
  df_b = 2^df_b
  colnames(df_b) = signatures_b[1]
}

df = cbind(df_a,df_b)



### build the integrated table
### single gene
colnames(df)[colnames(df) == signatures[1]] = symbol1
colnames(df)[colnames(df) == signatures[2]] = symbol2

pdf(file = outputdir,title="Result Display",width = 6,height = 5.5)
par(mar=c(4.5, 5.1, 1.1, 2.1))
options(warn=-1)
rvalue = signif(cor(x = df[,1], y = df[,2],method = method),2)
cpvalue = signif(as.numeric(cor.test(x = df[,1], y = df[,2],method = method)[3]),2)
plot(x = log2(df[,1] + 1),y = log2(df[,2] + 1),main = NULL,cex.lab = 1.5,
     xlab = paste("log2(",colnames(df)[1]," TPM)",sep = ""),
     ylab = paste("log2(",colnames(df)[2]," TPM)",sep = ""), pch = 19,cex=0.5)
range_y = range(log2(df[,2] + 1))
text(x = max(log2(df[,1] + 1)) * 1.03,y = max(log2(df[,2] + 1)) - (range_y[2] - range_y[1])/15,labels = paste("p-value = ",as.character(cpvalue),"\nR = ",as.character(rvalue),sep=""),cex = 1.3,col = "black",pos = 2)
a = dev.off()