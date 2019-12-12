datadir <- "/home/danielred/Programs/endoevol/OUT/TRBEC3_1/"

library(ggplot2)
getPar <- function(parfile, parname){
  strsplit(grep(parname, readLines(  parfile ), value=T), " ")[[1]][2]
}
setwd(datadir)
#death = getPar("parameters.txt", "par_death")
data <- read.table("output_tail.txt", header=F, sep="\t", col.names= c("ttime", "cell", "length", "role", "seq", "krepl", "kendo", "kT", "kC"), colClasses = c("character", NA, NA, "factor", "character", NA, NA, NA, NA))
data$ttime <- as.numeric(substring(data$ttime, 3))
#data <- cbind(data, length=nchar(data$seq))
str(data)
atlagadat <- data.frame(
time=rep(unique(data$ttime), 4),
number=c( sapply(unique(data$ttime), function(x, data) length(data[data$ttime == x, "krepl"]), data ),
         sapply(unique(data$ttime), function(x, data) length(data[data$ttime == x, "kendo"]), data ),
         sapply(unique(data$ttime), function(x, data) length(data[data$ttime == x, "kT"]), data ),
         sapply(unique(data$ttime), function(x, data) length(data[data$ttime == x, "kC"]), data )
),
activity=c( sapply(unique(data$ttime), function(x, data) mean(data[data$ttime == x, "krepl"]), data ),
            sapply(unique(data$ttime), function(x, data) mean(data[data$ttime == x, "kendo"]), data ),
            sapply(unique(data$ttime), function(x, data) mean(data[data$ttime == x, "kT"]), data ),
            sapply(unique(data$ttime), function(x, data) mean(data[data$ttime == x, "kC"]), data )
        ),
sd=c( sapply(unique(data$ttime), function(x, data) sd(data[data$ttime == x, "krepl"]), data ),
            sapply(unique(data$ttime), function(x, data) sd(data[data$ttime == x, "kendo"]), data ),
            sapply(unique(data$ttime), function(x, data) sd(data[data$ttime == x, "kT"]), data ),
            sapply(unique(data$ttime), function(x, data) sd(data[data$ttime == x, "kC"]), data )
        ),
type=rep(c("krepl", "kendo", "kT", "kC"), each=length(unique(data$ttime)) )
)
atlagadat <- cbind(atlagadat, se=atlagadat$sd/sqrt(atlagadat$number))

roletime <- data.frame(
  time= rep(unique(data$ttime), 5),
  number= c(
    sapply(unique(data$ttime), function (x, data) nrow(data[data$ttime==x & data$role %in% c(1, "single"),]), data),
    sapply(unique(data$ttime), function (x, data) nrow(data[data$ttime==x & data$role %in% c(2, "repl"),]), data),
    sapply(unique(data$ttime), function (x, data) nrow(data[data$ttime==x & data$role %in% c(3, "endo"),]), data),
    sapply(unique(data$ttime), function (x, data) nrow(data[data$ttime==x & data$role %in% c(4, "repl_template"),]), data),
    sapply(unique(data$ttime), function (x, data) nrow(data[data$ttime==x & data$role %in% c(5, "endo_template"),]), data)
  ),
  type= rep(c("single", "repl", "endo", "repl_template", "endo_template"), each=length(unique(data$ttime)) ),
  length= c(sapply(c("single", "repl", "endo", "repl_template", "endo_template"), function(x) sapply(unique(data$ttime), function (y, data) mean(data[data$ttime==y & data$role == x,"length"]), data)))
)


## FIGURES

table(data$role)
ggplot(atlagadat, aes(x=time, y=activity, colour=type))+
  geom_line()+
  geom_errorbar(aes(ymin=activity-se, ymax=activity+se))+
  labs(caption=paste("par_death", getPar("parameters.txt", "par_death"), "diff", getPar("parameters.txt", "par_diffusion_rate"), "par_substitution", getPar("parameters.txt", "par_substitution") ) )

ggplot(roletime, aes(x=time, y=number, fill=type))+
  geom_hline(yintercept=512*512)+
  geom_hline(yintercept=512*512/2)+
  geom_area()+
  labs(caption=paste("par_death", getPar("parameters.txt", "par_death"), "diff", getPar("parameters.txt", "par_diffusion_rate"), "par_substitution", getPar("parameters.txt", "par_substitution") ), title=getPar("parameters.txt", "par_ID") )
  
ggplot(roletime, aes(x=time, y=length, color=type))+
  geom_line()+
  labs(caption=paste("par_death", getPar("parameters.txt", "par_death"), "diff", getPar("parameters.txt", "par_diffusion_rate"), "par_substitution", getPar("parameters.txt", "par_substitution") ), title=getPar("parameters.txt", "par_ID") )

plot(data$ttime, data$length)
  
max(apply(data[,c("krepl","kendo", "kT", "kC")], 1, sum))
max(data[,c("krepl","kendo", "kT", "kC")])

hist(data[data$ttime==max(data$ttime),"kC"])
hist(data[data$ttime==max(data$ttime),"kT"])
hist(data[data$ttime==max(data$ttime),"krepl"])
hist(data[data$ttime==max(data$ttime),"length"])

    
