bmrROC <- function(ics,posVisit,negVisit,polyCytokine,stimulation,control,Parent,column="visit",margin=2,negate=FALSE,mcmc=FALSE,subset.mapping=NULL)
{
	filter<-substitute(filter)
	pos<-extractDataPoly(ics=ics,cytokineA="foo",cytokineB="foo",and=polyCytokine,stim=stimulation,control=control,subset=c(posVisit,Parent),subset.mapping=subset.mapping)
	neg<-extractDataPoly(ics=ics,cytokineA="foo",cytokineB="foo",and=polyCytokine,stim=stimulation,control=control,subset=c(negVisit,Parent),subset.mapping=subset.mapping)
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
	filter<-filterTcells(m)&filterBackground(m)
	m<-m[filter,]
	m<-merge(m,rx,"ID")
	mm<-m[,c("ns","Ns","nu","Nu")]
	attr(mm,"cytokine")<-c(polyCytokine)
	attr(mm,"class")<-c("data.frame")
	attr(mm,"control")<-c(control)
	attr(mm,"stimulation")<-c(stimulation)
	
	bmr<-BetaMix(mm,alternative="greater",scl=10)
	pData(bmr)<-m[,c("ID","rx_code","visit")]
	rocdata<-lapply(seq(0,1,l=500),function(i){
				d<-fisherVsBB(bmr,threshold=i,column="visit",subset=!rx_code%in%c("P1-P2","P3-P4"),margin=margin,adjust.fisher="none",adjust.bb="none")
				rn<-rownames(d[[1]])
				data.frame(melt(d[[1]]),visit=rn,threshold=i)
			} )
	tmp<-do.call(rbind,rocdata)
	colnames(tmp)<-c("Method","Response.Rate","Visit","Nominal.FDR")  
	
	tmp$Visit<-factor(tmp$Visit)

	Title1<-paste(stimulation," / ",polyCytokine," CD4+ Tcells from \n (visit ",negVisit,") and (visit ",posVisit,")",sep="")
	form<-as.formula(paste("`",posVisit,"`~`",negVisit,"`",sep=""))
	
	#compute the bonferroni corrected p-value and response rate and AUC
	thresh<-0.05/length(setdiff(levels(ics@antigen),c("CMV","sebctrl")))
	bonf<-data.frame(response=MIMOSA:::fisherTest(bmr,threshold=thresh,adjust="none"),visit=get("visit",pData(bmr)))
	bonf$response<-factor(bonf$response,levels=c("FALSE","TRUE"))
	subs<-with(pData(bmr),!rx_code%in%c("P1-P2","P3-P4"))
	bonf<-prop.table(table(bonf[subs,]),margin=2)["TRUE",]
	recast<-cast(melt(id=c("Method","Visit","Nominal.FDR"),measure="Response.Rate",tmp),Nominal.FDR+Method~Visit)
	myarg=c("TPR","FPR")
	attr(myarg,"names")<-c(posVisit,negVisit)
	recast<-(rename(recast,myarg))
	bonf<-bonf[c(posVisit,negVisit)]
	foo<-rename(bonf,myarg)
	foo<-as.data.frame(t(foo))
	
	l<-paste("AUC =",round(with(subset(recast,Method=="BB"),trapz(x=FPR,y=TPR)),digits=3))
	l2<-paste("AUC =",round(with(subset(recast,Method=="Fisher"),trapz(x=FPR,y=TPR)),digits=3))
	
	plt<-ggplot()+geom_path(data=recast,aes(x=FPR,y=TPR,group=Method,col=Method))+
			geom_point(data=foo,aes(x=FPR,y=TPR,size=2),shape=2)+
			scale_size_continuous(name="Bonferroni\nadjusted p-value",labels=round(thresh,digits=5)) + 
			geom_text(x=0.45,y=0.5,label=l,aes(color="BB"))+
			geom_text(x=0.35,y=0.4,label=l2,aes(color="Fisher"))+
			theme_bw()+opts(title=Title1)
	
	plt2<-ggplot()+geom_path(data=subset(tmp,Visit==2),aes(x=Nominal.FDR,y=Response.Rate,group=Method,col=Method))+theme_bw()+geom_abline(lty=2)+opts(title=Title1)+geom_hline(data=foo,aes(yintercept=FPR,alpha=1,color="Fisher"),lty=4,show_guide=TRUE)+scale_alpha(name="Observed False\nPositive Rate",label="Fisher's test at \nbonferroni\nadjusted\np-value")
	
	if(mcmc){
		d<-list(n.stim=bmr@data[,c("Ns","ns")],n.unstim=bmr@data[,c("Nu","nu")])
		inits<-MDMix(d,initonly=TRUE)
		r<-.fitMCMC(d,inits,100000,50000,4,alternative="not equal")
		d2<-t(sapply(seq(1,0,l=500),function(th)prop.table(table(factor(r$z[,2]>=th,levels=c("FALSE","TRUE"))[subs],pData(bmr)$visit[subs]),margin=2)["TRUE",]))
		colnames(d2)<-c("TPR","FPR")
		d2<-data.frame(d2,Method="MCMC two-sided")
		plt<-plt+geom_line(data=data.frame(d2),aes(x=FPR,y=TPR,color=Method))+geom_text(aes(color="MCMC two-sided"),x=0.65,y=0.75,label=paste("AUC =",round(trapz(d2[,"FPR"],d2[,"TPR"]),digits=3)))
	}
	#Fit two-sided Beta-binomial
	bmr.ne<-BetaMix(bmr@data,alternative="not equal")
	d3<-t(sapply(seq(1,0,l=500),function(th) prop.table(table(factor(bmr.ne@z[,2]>=th,levels=c("FALSE","TRUE"))[subs],pData(bmr)$visit[subs]),margin=2)["TRUE",]))
	colnames(d3)<-c("TPR","FPR")
	d3<-data.frame(d3,Method="BB two-sided")
	plt<-plt+geom_line(data=d3,aes(x=FPR,y=TPR,color=Method))+geom_text(aes(color="BB two-sided"),x=0.8,y=0.3,label=paste("AUC =",round(trapz(d3[,"FPR"],d3[,"TPR"]),digits=3)))
	
	#one sided mcmc

	oneside<-apply(bmr@data[,c("Ns","ns")],1,prop.table)[2,]<apply(bmr@data[,c("Nu","nu")],1,prop.table)[2,]
	r<-.fitMCMC(d,inits,100000,50000,4,alternative="greater")
	d4<-t(sapply(seq(1,0,l=500),function(th)prop.table(table(factor(r$z[,2]>=th,levels=c("FALSE","TRUE"))[subs],pData(bmr)$visit[subs]),margin=2)["TRUE",]))
	colnames(d4)<-c("TPR","FPR")
	d4<-data.frame(d2,Method="MCMC one-sided")
	plt<-plt+geom_line(data=data.frame(d4),aes(x=FPR,y=TPR,color=Method))+geom_text(aes(color="MCMC one-sided"),x=0.65,y=0.75,label=paste("AUC =",round(trapz(d4[,"FPR"],d4[,"TPR"]),digits=3)))
	
	return(list(roc=plt,fdr=plt2,bmr=bmr))
}



extractDataPoly<-function(ics = NULL, cytokineA = NULL, cytokineB = NULL, or = NULL,and=NULL,
		stim = NULL, control = NULL, subset = NULL, shrink = 1, scl = 1,...){
	
	if (any(c(is.null(ics), is.null(stim), 
					is.null(and), is.null(control), is.null(subset)))) {
		stop("All arguments must be non--null")
	}
	if(!is.null(or)){
		OR <- extractData(ics, control = control, stim = stim, subset = c(subset, 
						or),...)
		stop("Not supported yet")
	}else if(!is.null(and)){
		AND <- extractData(ics, control = control, stim = stim, subset = c(subset, and),...)
		return(AND)
	}
	A <- extractData(ics, control = control, stim = stim, subset = c(subset, 
					cytokineA),...)
	B <- extractData(ics, control = control, stim = stim, subset = c(subset, 
					cytokineB),...)		
}
