data(ICS)

#debugonce(ConstructMIMOSAExpressionSet)

E<-ConstructMIMOSAExpressionSet(ICS,
   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
   .variables= .(TCELLSUBSET,CYTOKINE,UID),
   featureCols=1,ref.append.replace='_REF')




