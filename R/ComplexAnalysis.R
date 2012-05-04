bmrROC <- function(ics,posVisit,negVisit,polyCytokine,stimulation,control,Parent,column="visit")
{
	filter<-substitute(filter)
	pos<-extractDataPoly(ics=ics,cytokineA="foo",cytokineB="foo",and=polyCytokine,stim=stimulation,control=control,subset=c(posVisit,Parent))
	neg<-extractDataPoly(ics=ics,cytokineA="foo",cytokineB="foo",and=polyCytokine,stim=stimulation,control=control,subset=c(negVisit,Parent))
	pos$ID<-rownames(pos)
	neg$ID<-rownames(neg)
	
	neg$visit=negVisit
	pos$visit=posVisit
	m<-rbind(neg,pos)
	
	merged<-mergeICSData(ics)
	pd<-subset(merged,ID%in%rownames(m)&cytokine%in%polyCytokine&antigen%in%stimulation&(visit%in%negVisit|visit%in%posVisit)&parent%in%Parent)
	pd<-cast(melt(pd,id=c("ID","visit")))
	pd<-unique(subset(pd,select=c("ID","rx_code")))
	rx<-pd
	
	m<-merge(m,rx,"ID")
	mm<-m[,c("ns","Ns","nu","Nu")]
	attr(mm,"cytokine")<-c(polyCytokine)
	attr(mm,"class")<-c("data.frame")
	attr(mm,"control")<-c(control)
	attr(mm,"stimulation")<-c(stimulation)
	
	bmr<-BetaMix(mm,alternative="greater",scl=10)
	pData(bmr)<-m[,c("ID","rx_code","visit")]
	rocdata<-mclapply(seq(0,1,l=500),function(i){
				d<-fisherVsBB(bmr,threshold=i,column="visit",subset=!rx_code%in%c("P1-P2","P3-P4"))
				rn<-rownames(d[[1]])
				data.frame(melt(d[[1]]),visit=rn,threshold=i)
			} )
	tmp<-do.call(rbind,rocdata)
	colnames(tmp)<-c("Method","Response.Rate","Visit","Nominal.FDR")   
	tmp$Visit<-factor(tmp$Visit)
	#tmp$Visit<-factor(tmp$Visit,labels=c("Day 0", "50 Days Post Vaccine 3"))
	Title<-paste(stimulation," / ",polyCytokine," CD4+ Tcells from \n (visit ",negVisit," ) and (visit ",posVisit,")")
	form<-as.formula(paste("`",posVisit,"`~`",negVisit,"`",sep=""))
	#xyplot(`Response.Rate`~`Nominal.FDR`|Visit,group=Method,tmp,auto.key=T,type="l",lwd=2,main=T,par.strip.text=list(cex=0.8),par.settings=list(par.main.text=list(cex=0.8)))
	xyplot(`12`~`2`,group=Method,cast(melt(id=c("Method","Visit","Nominal.FDR"),measure="Response.Rate",tmp),Nominal.FDR+Method~Visit),main=Title,par.settings=list(par.main.text=list(cex=0.8)),type="l",auto.key=TRUE,lwd=2,xlab="FPR",ylab="TPR")
}



extractDataPoly<-function(ics = NULL, cytokineA = NULL, cytokineB = NULL, or = NULL,and=NULL,
		stim = NULL, control = NULL, subset = NULL, shrink = 1, scl = 1){
	
	if (any(c(is.null(ics), is.null(stim), 
					is.null(and), is.null(control), is.null(subset)))) {
		stop("All arguments must be non--null")
	}
	if(!is.null(or)){
		OR <- extractData(ics, control = control, stim = stim, subset = c(subset, 
						or))
		stop("Not supported yet")
	}else if(!is.null(and)){
		AND <- extractData(ics, control = control, stim = stim, subset = c(subset, and))
		return(AND)
	}
	A <- extractData(ics, control = control, stim = stim, subset = c(subset, 
					cytokineA))
	B <- extractData(ics, control = control, stim = stim, subset = c(subset, 
					cytokineB))		
}
