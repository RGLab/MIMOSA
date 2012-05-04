#extractData.R

# Created on: Nov 29, 2011
#     Author: Greg Finak
#S4 CLASSES



#Function that sets up the matrix to contrast a stimulation with a negative control from ICS data.
#The input is from the ICS pipeline and containst columns "pos","neg","antigen","fname","ID","fcsfile","MFI","parent","pub_id"
#parent is the parent population for the cytokine. 
#fcsfile is the originating fcs file name
#antigen is the stimulation used for the observation
#pos and neg are the positive and negative cell counts
#ID is the sample id
#MFI is the median fluorescence intensity (not used here)
#pub_id is another unique sample identifier
#fname is the cytokine name
#If intial.contrasts==TRUE, then restructure the matrix above to something more suitable using reshape package.
#The data will be turned into a nested list, nested by parent, fname, ID and contain a matrix of antigen, pos, neg for each ID
#If initial.cast==TRUE, we do this restructuring in the extractData function call.
#control is the name of the negative control antigen
#stim is the name of the positive stimulation antigen
#default values are provided for control and stim that may not work with all data
#subset are the levels of parent and fname to be extracted from the nesting. Order is important

#TODO update extractData to deal with double-positive cytokines.
setGeneric("extractData",function(ics,control,stim,subset){ standardGeneric("extractData")});
setMethod("extractData",c("ICS","character","character","character"),
		function(ics,control="Neg Cont",stim="CMV",subset=c("cd4","IFNg")){
			x<-ics@.Data
			e<-eval(parse(text=eval(paste('x[["',paste(subset,collapse='"]][["'),'"]]',sep=""))))
			if(is.null(e)){
				stop("The levels for nesting in subset are not valid: provided: ",paste(subset,collapse=","))
			}
			#organize the data frame so that we have posivie and negative counts for stimulation and control along the columns, then combine into one matrix for all subjects, given the cytokine and parent population combination
			#Only returns data if the subgroup has both the control and stimulated antigen
			myframe<-do.call(rbind,lapply(e,function(x){
								if(!is.null(x)){
								sub<-subset(x,antigen==stim|antigen==control)
								if(nrow(sub)!=2){
									return(NULL)
								}
								cast(melt(as.data.frame(sub),id="antigen"),~antigen+variable)
							}
							}))
						
			#Drop the "value" column, we don't need it. 
			myframe<-subset(myframe,select=setdiff(colnames(myframe),"value"))
			#construct new column names that are required by MIMOSA code "ns","nu", "Ns","Nu", the positive (n) and negative (N) stimulated (s) and unstimulated (u) counts.
			oldnames<-apply(expand.grid(c(stim,control),c("pos","neg")),1,function(x)paste(x,collapse="_"))
			myframe<-myframe[,oldnames]
			colnames(myframe)<-c("ns","nu","Ns","Nu")
			attr(myframe,"stimulation")<-stim
			attr(myframe,"control")<-control
			attr(myframe,"cytokine")<-subset
			#TODO construct an AnnotatedDataFrame for each observation
			return(myframe)
		})


