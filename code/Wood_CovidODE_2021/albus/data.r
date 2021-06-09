## The data appear to be provided in data_rtm.csv + a serology file.
## most of the columns of data_rtm.csv have been set to NA

# death3 is hospital deaths
# death 2 - death 3 is carehome deaths
# phe_admissions is hospital admissions
# phe_occupied_mv_beds is ICU ocupancy
# phe_patients is hospital bed occupancy.
# react_positive and react_samples are react PCR survey data.
#library(sircovid)
#for (loc in region) { ch <- carehomes_parameters(1,loc)
# fn <- paste("ch-",loc,".rda",sep="")
# save(ch,file=fn)
#}


## Sloppily this picks up location from environment
## called from...

cat(location,"\n") 
load(paste("ch-",location,".rda",sep=""))

library(timeDate)
## example of use... when was lockdown
julian(timeDate("2020-03-24"),origin=timeDate("2019-12-31"))

rtm <- read.csv("data_rtm.csv")
early <- read.table("early.dat",header=TRUE)
early$day <- as.integer(julian(timeDate(early$date,format="%d/%m/%Y"),
             origin=timeDate("2019-12-31")))
early <- early[,-(2:7)] ## drop national
early$north_east_and_yorkshire <- early$NE+early$York
early$midlands <- early$EMid+early$Wmid
early <- early[,-c(2,4,5,6)]
names(early) <- c("date","north_west","east_of_england","london","south_east","south_west",
"day","north_east_and_yorkshire","midlands")

sero <- read.csv("data_serology.csv")
sero$day <- as.integer(julian(timeDate(sero$date),origin=timeDate("2019-12-31")))
colSums(!is.na(rtm)) -> cs
rtm <- rtm[,cs>0] ## retain only 12/39 cols containing something
region <- unique(rtm$region) ## 7 regions


ser <- sero[sero$region==location,]
rtr <- rtm[rtm$region==location,]
rtr$day <- 75+1:nrow(rtr)
rtr$n_positive <- rtr$total_samples <- NA
ii <- match(ser$day,rtr$day)
rtr$n_positive[ii] <- ser$n_positive
rtr$total_samples[ii] <- ser$total_samples

## append earlier death records
if (FALSE) {
  dstart <- 50
  day <- dstart:75
  death <- day*0
  ii <- early$day>=dstart
  eday <- early$day[ii];edeath <- early[,location][ii]
  ii <- match(eday,day)
  death[ii] <- edeath
  stub <- as.data.frame(matrix(NA,length(death),ncol(rtr)))
  names(stub) <- names(rtr)
  stub$day <- day;stub$death3 <- death
  rtr <- rbind(stub,rtr)
}

## the rep41 b(t) knots
bk <- c("2020-03-16","2020-03-23","2020-03-25","2020-05-11","2020-06-15","2020-07-04",
   "2020-08-03","2020-09-01","2020-09-14","2020-10-14","2020-10-31","2020-11-05")
bk <- as.integer(julian(timeDate(bk),origin=timeDate("2019-12-31")))




