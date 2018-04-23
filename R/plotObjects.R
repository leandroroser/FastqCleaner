
#' plotObjects
#' Create the information required to construct the plots. 
#' This is the input of myplot, which uses the values created for 
#' this function to construct the plots
#' @return List with information to construct the diagnostic plots
#' @keywords internal

plotObjects <- function(fq, klength, basename, maxFreq, sampleSize) { 

out <- list()
fq.width <- width(fq)
fq.length <- length(fq)


if(fq.length == 0) {
    return(out)
}

input_seqs <- sread(fq)
seqs <- input_seqs
qual <- PhredQuality(quality(fq))


if(length(unique(fq.width)) == 1) {
    out$qual <- as.matrix(qual) 
    out$seq <- as.matrix(seqs)
} else {


mymax <- max(fq.width)
cuales <- which(width(seqs) < mymax)
agregar <- mymax - width(seqs[cuales])


to_paste_s <- cpp_create_stringvec(agregar, "-")
seqs[cuales] <- paste0(seqs[cuales], to_paste_s)

to_paste_q <- cpp_create_stringvec(agregar, " ")
qual[cuales] <- paste0(qual[cuales], to_paste_q)

out$seq <-  as.matrix(seqs)
out$qual <- as.matrix(qual)

out$seq[out$seq=="-"] <- NA
mode(out$qual) <- "numeric"
out$qual[out$qual == -1] <- NA
}
colnames(out$qual) <- colnames(out$seq) <- seq_len(ncol(out$seq))




bprange <- apply(out$qual, 2, function(x) { 
round(quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9),  na.rm = TRUE))
})
out$bp_cycle <- data.frame(colnames(out$qual), t(bprange))
colnames(out$bp_cycle) <- c("Cycle", "min", "low", "mid", "top", "max")


A_ <- T_ <- G_ <- C_ <- out$qual

A_[out$seq %in% c("T", "G", "C")] <- NA
A_ <- colMeans(A_, na.rm=TRUE)

T_[out$seq %in% c("A", "G", "C")] <- NA
T_ <- colMeans(T_, na.rm=TRUE)

G_[out$seq %in% c("T", "A", "C")] <- NA
G_ <- colMeans(G_, na.rm=TRUE)

C_[out$seq %in% c("T", "G", "A")] <- NA
C_ <- colMeans(C_, na.rm=TRUE)

out$meanBaseQuality <- data.frame(Quality=c(A_, C_, G_, T_), 
Base=rep(c("A", "C", "G", "T"), 
        each=length(A_)), 
        Cycle=rep(colnames(out$qual, 4)))

out$meanBaseQuality[, 1] <- round(out$meanBaseQuality[, 1], 2)



qv <- table(round(rowMeans(out$qual)))[as.character(seq(0, max(out$qual,
na.rm=TRUE), 1))]
qv[is.na(qv)] <- 0
names(qv) <- seq(0, max(out$qual, na.rm=TRUE), 1)
meanQReads <- data.frame(qv)
colnames(meanQReads) <- c("Quality", "Percent")
meanQReads[,2] <- 100 * meanQReads[,2] / sum(meanQReads[, 2])
meanQReads[,1] <- factor(meanQReads[,1], 
                    levels=unique(meanQReads[,1]), 
                    ordered=TRUE)
                    out$meanQReads <- meanQReads


                    
qualth <- seq_len(40)
qual_above <- cpp_check_quality_threshold(out$qual, qualth)
qual_above <- 100 * colSums(qual_above) / nrow(out$qual)
out$qual_above <- qual_above

base_content <- 100 * cpp_base_content(out$seq)
out$base_content <- round( base_content , 2)

GC_Content <- cpp_GC_content(out$seq)
GC_Content <- 100 * rowSums(GC_Content[,seq_len(2)],
                        na.rm = TRUE) / rowSums(GC_Content, na.rm =  TRUE)
out$GC <- hist(GC_Content, plot = FALSE, breaks = "Sturges")



freqreads <- table(fq.width)
out$length_distribution$ydata <- unname(freqreads / sum(freqreads) * 100)
out$length_distribution$xdata <- factor(names(freqreads), 
levels=unique(names(freqreads)), ordered=TRUE)



qa1 <- ShortRead::qa(fq, basename(basename), n = sampleSize)
seqtable <- qa1[["frequentSequences"]]
counts <- as.numeric(as.vector(seqtable$count)) 
myDT <-data.frame(round(100* counts/sampleSize, 1))
myDT <-data.frame(Sequences = seqtable[, 1], 
Percent = round(100*seqtable[, 2]/sampleSize, 1))
colnames(myDT)[2] <- paste0("% of sample (n=", sampleSize, ")")
out$dt <- myDT

myintervals <- data.frame(labels = c('1', '2-10', '11-100', 
'101-1k', '1k-10k', '>10k'), 
lower = c(1,2,11,101,1001,10001), 
upper = c(2,11,101,1001,10001, Inf))

count_ocurrence <- vapply(seq_along(myintervals[,1]), 
                function(x) {
                sum(counts >= myintervals[x,2] & counts < myintervals[x,3])
                }, integer(1))
ocurrenceData <- data.frame(labels = myintervals[,1], 
Percent = count_ocurrence)
ocurrenceData[, 1] <- factor(ocurrenceData[,1], 
                        levels = unique(ocurrenceData[,1]), 
                        ordered = TRUE)
                        out$ocurrenceData <- ocurrenceData


                        
                        
oligofreq <- oligonucleotideFrequency(input_seqs, klength)
oligofreq <- colSums(oligofreq)
oligosum <- sum(oligofreq)
oligofreq <- oligofreq[oligofreq  != 0]
oligofreq <- oligofreq[order(oligofreq , decreasing = TRUE)]
oligofreq <- round(100 * oligofreq/oligosum, 1)
oligofreq <- data.frame(sequence = names(oligofreq), 
percent = unname(oligofreq))
colnames(oligofreq)[2]<-"% in sample"
out$oligofreq <- oligofreq
out$oligosum <- oligosum
out$klength <- klength


loopv <- seq_len(min(width(seqs) - (klength-1)))

kcount <- lapply(loopv, function(y) { 
DNAStringSet(start=y, end=y+klength-1, seqs)
})

if(length(unique(fq.width)) == 1) {
    kcount <- vapply(kcount, function(y) length(unique(y)), integer(1))
} else {


kcount <- vapply(kcount, function(y) {
cuales <- grep("-", y)
if(length(cuales) == 0) 
    return(length(unique(y)))
else
    return(length(unique(y[-grep("-", y)])))
}, integer(1))
}
reldiv <- kcount/(5^klength) 
reldiv <- c(reldiv, rep(NA, klength-1)) 
names(reldiv) <- seq_along(reldiv)
kmerData <- data.frame(RelDiv = reldiv, 
Method = rep(c(1), 
each = length(reldiv)), 
Cycle = names(reldiv))
out$kmerData  <- kmerData


out

}


#' myPlot
#' @description Construction of diagnostic plots. The function depends
#' of the values created by plotObject
#' @return List with Highcharts plots
#' @keywords internal


myPlot <- function(isPaired, location, sampleSize,
kmerLength, theFile, maxFreq) {

plotList <- function(out, nplots, theFile, sampleSize) {

if(length(out) == 0) {
    out <- alist(none=, a=, b=, c=, d=, e=, f=,
            g=, h=, i=, j=,table=, oligofreq=)
for(i in seq_len(11)) {
    out[[i]] = shiny::HTML("<script></script>")
}
out[[12]] <- out[[13]] <-  data.frame()
return(out)
}

list(none = "<script></script>",
a = plotA(out, nplots, theFile, sampleSize),
b = plotB(out, nplots, theFile, sampleSize), 
c = plotC(out, nplots, theFile, sampleSize), 
d = plotD(out, nplots, theFile, sampleSize), 
e = plotE(out, nplots, theFile, sampleSize),
f = plotF(out, nplots, theFile, sampleSize), 
g = plotG(out, nplots, theFile, sampleSize),
h = plotH(out, nplots, theFile, sampleSize),
i = plotI(out, nplots, theFile, sampleSize),
j = plotJ(out, nplots, theFile, sampleSize),
table = out$dt, 
oligofreq = out$oligofreq)
}  


if(!isPaired) {
stream <- FastqSampler(location, sampleSize)
data <- yield(stream)
x <- list(plotList(plotObjects(data,  kmerLength, location,
                maxFreq, sampleSize), 
                nplots = 1, theFile, sampleSize))
close(stream)
} else {
streamL <- FastqSampler(location[1], sampleSize)
dataL <- yield(streamL)
streamR <- FastqSampler(location[2], sampleSize)
dataR <- yield(streamR)
x <- list(plotList(plotObjects(dataL, kmerLength, location[1], 
                    maxFreq, sampleSize),  nplots = 1, theFile, sampleSize),
                    plotList(plotObjects(dataR, kmerLength, location[2], 
maxFreq, sampleSize),  nplots = 2, theFile, sampleSize))
close(streamL)
close(streamR)
}
x
}


#' plotA
#' @return Per cycle quality  plot
#' @keywords internal


plotA <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
container <- "container1"
} else  {
container <- "container2"
}

theFile <- match.arg(theFile)

jsondata <- lapply(unname(split(x$bp_cycle[,-1],seq_len(nrow(x$bp_cycle)))),
                    function(y)paste(y, collapse = ", "))
jsondata <- paste("[", jsondata, "]", sep = "", collapse = ",\n")
jsondata <- paste0("[", jsondata, "]")

xAxis <- paste("[",  paste("'", x$bp_cycle[, 1], "'", sep = "", 
collapse = ", "), "]", sep ="")

ymax <- max(x$bp_cycle[,-1], na.rm = TRUE)

if(ymax < 40) { 
    ymax <- 40
}

ymin <- min(x$bp_cycle[,-1], na.rm = TRUE)

if(ymin > 0) { 
    ymin <- 0
}


paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({

chart: {
type: 'boxplot'
},

title: {
text: 'Per cycle quality - ", theFile, "'
},

legend: {
enabled: false
},

xAxis: { 
categories:", xAxis,",
title: {
text: 'Cycle'
}
},

yAxis: {
max:", ymax, ",
min:", ymin, ",
title: {
text: 'Observations'
}
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
plotOptions: {
boxplot: {
fillColor: '#95CEFF',
lineWidth: 0,
medianColor: '#2A2A93',
medianWidth: 3,
stemColor: '#A63400',
stemDashStyle: 'dot',
stemWidth: 1,
whiskerColor: '#3D9200',
whiskerLength: '50%',
whiskerWidth: 3
}
},

series: [{
name: 'Observations',
data:", jsondata,"}],
exporting: {
enabled: true
},

credits: {
enabled: false
}

});
});
</script>")
}


#' plotB
#' @return Per cycle mean base quality plot
#' @keywords internal


plotB <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {

if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


xAxis <- paste("[",  paste("'", x$meanBaseQuality[, 3], "'", sep = "", 
collapse = ", "), "]", sep ="")
ymax <- round(max(x$meanBaseQuality[,1], na.rm = TRUE))
if(ymax < 40) { 
    ymax <- ymax  + 1
}

meanBaseQuality <- split(x$meanBaseQuality[, 1], x$meanBaseQuality[,2])
nombres <- names(meanBaseQuality)
jsondata <- lapply(unname(meanBaseQuality),  
                    function(x)paste(x, collapse = ", "))
jsondata <- paste("[", jsondata, "]", sep = "")
names(jsondata) <- nombres


paste0("<script type= 'text/javascript'>

$(function () {
$('#", container, "').highcharts({
title: {
x: -20, //center
text: 'Per cycle mean base quality - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
"},
yAxis: {
max:", ymax, ",
title: {
text: 'Quality'
},
plotLines: [{
value: 0,
width: 1,
color: '#808080'
}]
},
tooltip: {
pointFormat: '<span style=\"color:{series.color}\">{series.name}</span>:'+ 
'<b>{point.y}%</b><br/>'
},
legend: {
layout: 'vertical',
align: 'right',
verticalAlign: 'middle',
borderWidth: 0
},
credits: {
enabled: false
},
series: [{
name:",  paste0("\"",names(jsondata)[1],"\""), ",\n",
"data:", jsondata[[1]],"\n",
"}, {
name:",   paste0("\"",names(jsondata)[2],"\""), ",\n",
"data:", jsondata[[2]],"\n",
"}, {
name:",   paste0("\"",names(jsondata)[3],"\""), ",\n",
"data:", jsondata[[3]], "\n",
"}, {
name:",   paste0("\"",names(jsondata)[4],"\""), ",\n",
"data:", jsondata[[4]], "\n",
"}]
});
});
</script>")

}


#' plotC
#' @return Mean quality of reads distribution plot
#' @keywords internal

plotC <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)

x$meanQReads[, 2] <- round(x$meanQReads[, 2], 2) 
yrange <- range(x$meanQReads[, 2], na.rm = TRUE)
ymax <- min(yrange[2] + yrange/10, 100)

xAxis <- paste("[", paste("'", x$meanQReads[, 1], "'", 
collapse = ", ", sep = ""), "]", sep ="")

plotReads <- unname(x$meanQReads[, 2])
jsondata <- paste(plotReads, collapse = ", ")
jsondata <- paste("[", jsondata, "]", sep = "")


paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({
chart: {
type: 'column'
},
title: {
text: 'Distribution of mean quality of reads - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
",
crosshair: true,
title: {
text: 'Mean Quality'
}
},
yAxis: {
min: 0,
max:", ymax, ",
title: {
text: '% Reads'
}
},
tooltip: {
pointFormat: '<b>{point.y}%</b><br/>'
},
plotOptions: {
column: {
pointPadding: 0.2,
borderWidth: 0
}
},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
data:", jsondata,"\n",
"}]
});
});
</script>")

}


#' plotD
#' @return  percent of reads with quality > threshold plot
#' @keywords internal

plotD <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


jsondata <- paste("[", paste(x$qual_above, collapse = ", "), "]", sep = "")
xAxis <- paste("[",  paste("'", seq_len(40), "'", sep = "",
collapse = ", "), "]", sep ="")

ymin <- 0
ymax <- 100


paste0("<script type= 'text/javascript'>

$(function () {
$('#", container, "').highcharts({
title: {
x: -20, //center
text: 'Quality - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
"},
yAxis: {
max:", ymax, ",
min:", ymin, ",
title: {
text: '% reads'
},
plotLines: [{
value: 0,
width: 1,
color: '#808080'
}]
},
tooltip: {
pointFormat: '<b>{point.y}%</b><br/>'
},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
name:",  paste0("\"",names(jsondata),"\""), ",\n",
"data:", jsondata,"\n",
"}]
});
});
</script>")

}

#'plotE
#' @return Per cycle base proportion plot
#' @keywords internal


plotE <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


ymax <- 100
ymin <- 0

jsondata <- list(
'C' = paste("[", paste(x$base_content[1, ], 
collapse = ", "), "]", sep = ""),
'G' = paste("[", paste(x$base_content[2, ], 
collapse = ", "), "]", sep = ""),
'A' = paste("[", paste(x$base_content[3, ], 
collapse = ", "), "]", sep = ""),
'T' = paste("[", paste(x$base_content[4, ], 
collapse = ", "), "]", sep = ""),
'N' = paste("[", paste(x$base_content[5, ], 
collapse = ", "), "]", sep = "")
)


xAxis <- paste("[",  paste("'", names(x$base_content), "'", 
sep = "", collapse = ", "), "]", sep ="")


paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({
chart: {
type: 'column'
},
title: {
text: 'Per cycle base proportion (%) - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:",  xAxis, 
"},
yAxis: {
min: 0,
max: 100,
title: {
text: 'Base proportion (%)'
}
},
tooltip: {
pointFormat: '<span style=\"color:{series.color}\">{series.name}</span>:'+ 
'<b>{point.y}%</b><br/>',
shared: true
},
plotOptions: {
column: {
stacking: 'percent'
}
},
credits: {
enabled: false
},
series: [{
name:",  paste0("\"",names(jsondata)[1],"\""), ",\n",
"data:", jsondata[[1]],"\n",
"}, {
name:",   paste0("\"",names(jsondata)[2],"\""), ",\n",
"data:", jsondata[[2]],"\n",
"}, {
name:",   paste0("\"",names(jsondata)[3],"\""), ",\n",
"data:", jsondata[[3]], "\n",
"}, {
name:",   paste0("\"",names(jsondata)[4],"\""), ",\n",
"data:", jsondata[[4]], "\n",
"}, {
name:",   paste0("\"",names(jsondata)[5],"\""), ",\n",
"data:", jsondata[[5]], "\n",
"}]
});
});
</script>"
)

}


#'plotF
#'@return Per cycle base proportion plot (lineplot)
#' @keywords internal


plotF <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)

yrange <- range(x$base_content, na.rm = TRUE)
ymax <- min(yrange[2] + yrange[2]/10, 100)
ymin <- max(yrange[1] - yrange[1]/10, 0)

jsondata <- list(
'C' = paste("[", paste(x$base_content[1, ],
collapse = ", "), "]", sep = ""),
'G' = paste("[", paste(x$base_content[2, ], 
collapse = ", "), "]", sep = ""),
'A' = paste("[", paste(x$base_content[3, ], 
collapse = ", "), "]", sep = ""),
'T' = paste("[", paste(x$base_content[4, ], 
collapse = ", "), "]", sep = ""),
'N' = paste("[", paste(x$base_content[5, ], 
collapse = ", "), "]", sep = "")
)


xAxis <- paste("[",  paste("'", names(x$base_content), "'", sep = "", 
collapse = ", "), "]", sep ="")


paste0("<script type= 'text/javascript'>

$(function () {
$('#", container, "').highcharts({
title: {
x: -20, //center
text: 'Per cycle base proportion (%) - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
"},
yAxis: {
max:", ymax, ",
min:", ymin, ",
title: {
text: 'Base proportion (%)'
},
plotLines: [{
value: 0,
width: 1,
color: '#808080'
}]
},
tooltip: {
pointFormat: '<span style=\"color:{series.color}\">{series.name}</span>:'+ 
'<b>{point.y}%</b><br/>'
},
legend: {
layout: 'vertical',
align: 'right',
verticalAlign: 'middle',
borderWidth: 0
},
credits: {
enabled: false
},
series: [
{
name:",  paste0("\"",names(jsondata)[1],"\""), ",\n",
"data:", jsondata[[1]],"\n",
"}, 
{
name:",   paste0("\"",names(jsondata)[2],"\""), ",\n",
"data:", jsondata[[2]],"\n",
"},
{
name:",   paste0("\"",names(jsondata)[3],"\""), ",\n",
"data:", jsondata[[3]],"\n",
"},

{
name:",   paste0("\"",names(jsondata)[4],"\""), ",\n",
"data:", jsondata[[4]],"\n",
"},
{
name:",   paste0("\"",names(jsondata)[5],"\""), ",\n",
"data:", jsondata[[5]],"\n",

"}]
});
});
</script>")


}


#' plotG
#' @return CG content distribution plot
#' @keywords internal


plotG <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


ydata <- 100 * round(x$GC$counts/sum(x$GC$counts), 2)

yrange <- range(ydata)
ymax <- min(yrange[2] + yrange/10, 100)

jsondata <- paste("[", paste(ydata, collapse = ", "), "]", sep = "")
xAxis <- paste("[",paste( "'", x$GC$breaks,"'", 
sep = "", collapse = ", "), "]", sep ="")


paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({
chart: {
type: 'column'
},
title: {
text: '% reads with a given % GC - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
",
crosshair: true,
title: {
text: '% GC'
}
},
yAxis: {
min: 0,
max:", ymax, ",
title: {
text: '% Reads'
}
},
tooltip: {
pointFormat: '<b>{point.y}%</b><br/>'
},
plotOptions: {
column: {
pointPadding: 0.2,
borderWidth: 0
}
},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
data:", jsondata,"\n",
"}]
});
});
</script>")

}


#' plotH
#' @return Read length distribution
#' @keywords internal

plotH <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


yrange <- range(x$length_distribution$ydata, na.rm = TRUE)
ymax <- min(yrange[2] + yrange/10, 100)

xAxis <- paste("[", paste("'", x$length_distribution$xdata, "'", 
collapse = ", ", sep = ""), "]", sep ="")

jsondata <- paste(x$length_distribution$ydata, collapse = ", ")
jsondata <- paste("[", jsondata, "]", sep = "")


paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({
chart: {
type: 'column'
},
title: {
text: 'Read length distribution - ", theFile, "' 
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
",
crosshair: true,
title: {
text: 'Length (bp)'
}
},
yAxis: {
min: 0,
max:", ymax, ",
title: {
text: '% Reads'
}
},
tooltip: {
headerFormat: '<span style=\"font-size:10px\">{point.key}</span><table>',
pointFormat: '<tr><td style=\"color:{series.color};padding:0\">{series.name}:'+
'</td> <td style=\"padding:0\"><b>{point.y:.1f} % </b></td></tr>',
footerFormat: '</table>',
shared: true,
useHTML: true
},
plotOptions: {
column: {
pointPadding: 0.2,
borderWidth: 0
}
},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
data:", jsondata,"\n",
"}]
});
});
</script>")

}


#' plotI
#' @return Read ocurrence distribution plot
#' @keywords internal


plotI <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}
theFile <- match.arg(theFile)


yrange <- range(x$ocurrenceData[, 2], na.rm = TRUE)
ymax <- min(yrange[2] + yrange/10, 100)

xAxis <- paste("[", paste("'", x$ocurrenceData[, 1], "'", 
collapse = ", ", sep = ""), "]", sep ="")

ocurrenceData <- unname(x$ocurrenceData[,2])
jsondata <- paste(ocurrenceData, collapse = ", ")
jsondata <- paste("[", jsondata, "]", sep = "")



paste0("<script type= 'text/javascript'>
$(function () {
$('#", container, "').highcharts({
chart: {
type: 'column'
},
title: {
text: 'Read occurrence distribution - ", theFile, "' 
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
categories:", xAxis,
",
crosshair: true,
title: {
text: 'Read Ocurrence'
}
},
yAxis: {
min: 0,
max:", ymax, ",
title: {
text: '% Reads'
}
},
tooltip: {
headerFormat: '<span style=\"font-size:10px\">{point.key}</span><table>',
pointFormat: '<tr><td style=\"color:{series.color};padding:0\">{series.name}:'+
'</td> <td style=\"padding:0\"><b>{point.y:.1f} % </b></td></tr>',
footerFormat: '</table>',
shared: true,
useHTML: true
},
plotOptions: {
column: {
pointPadding: 0.2,
borderWidth: 0
}
},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
data:", jsondata,"\n",
"}]
});
});
</script>")

}



#' plotJ
#' @return Relative kmer diversity plot
#' @keywords internal

plotJ <- function(x, nplots = 1, theFile = c("input", "output"), sampleSize) {
if(nplots == 1) {
    container <- "container1"
} else  {
    container <- "container2"
}

theFile <- match.arg(theFile)


xAxis <- paste("[",  paste("'", x$kmerData[, 3], "'", sep = "", 
collapse = ", "), "]", sep ="")
kmerData <- split(x$kmerData[, 1], x$kmerData[, 2])
kmerData[[1]][is.na(kmerData[[1]])] <- NaN
nombres <- names(kmerData)
jsondata <- lapply(unname(kmerData), function(y)paste(y, collapse = ", "))
jsondata <- paste("[", jsondata, "]", sep = "")

names(jsondata) <- ""
ylab <- paste("k", x$klength, "-mer diversity", sep="")

paste0("<script type= 'text/javascript'>

$(function () {
$('#", container, "').highcharts({
title: {
x: -20, //center
text: 'Relative k-mer diversity (observed/expected k-mers) - ", theFile, "'
},
subtitle: {
text: 'sample size: ", sampleSize, "',
x: -20
},
xAxis: {
title: {
text:'Cycle' 
},
categories:", xAxis,
"},
yAxis: {
title: {
text:'", ylab,"' 
},
plotLines: [{
value: 0,
width: 1,
color: '#808080'
}]
},
//tooltip: {
//  valueSuffix: 'kmer-diversity'
//},
legend: {
enabled: false,
},
credits: {
enabled: false
},
series: [{
name:",  paste0("\'",names(jsondata),"\'"), ",\n",
"data:", jsondata,"\n",
"}]
});
});
</script>")

}
