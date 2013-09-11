[![Build Status](https://travis-ci.org/RGLab/MIMOSA.png?branch=master)](https://travis-ci.org/RGLab/MIMOSA)

# MIMOSA: Mixture Models for Single Cell Assays

MIMOSA is a package for fitting mixtures of beta-binomial or dirichlet-multinomial models to paired count data from single-cell assays, as typically appear in immunological studies (i.e. ICS, intracellular cytokine staining assay, or Fluidigm Biomark single-cell gene expression assays). 

The method is, generally, more sensitive and specific to detect differences between conditions (i.e. stimulated vs. unstimulated samples) than alternative approaches such as Fisher's exact test, or empirical ad-hoc methods like ranking by log-fold change.

I've cleaned up the package significantly. Old code has been rewritten, hidden, and the relevant stuff has been documented. The user-facing interface has been unified to use the 'MIMOSA' function to fit data. Data is represented by 'ExpressionSet' objects from 'Biobase'. See the vignette for an example.

MIMOSA shares information across subjects by means of priors on the proportions of stimulated and unstimulated cells, respectively.

### NEWS
- 09/10/2013 - MIMOSA() now takes default values ref=RefTreat%in%"Reference", subset=RefTreat%in%"Treatment". You no longer need to specify it explicitly. If it is missing from your own call, MIMOSA will warn you and add the additional factors. To disable the warning, set RT=FALSE.
- 09/06/2013 - Empty combinations of levels of conditioning variables are now dropped.
- Support for mclapply via parallel package in favor of older multicore package. Works with R 3.0.0 and 3.0.1
- 05/08/2013 - Compilation bugs fixed for Windows systems. 

### INSTALLATION

```r
install.packages("devtools") #install devtools package
library(devtools) #load it
install_github("MIMOSA","RGLab",branch="master") #Need R 3.0.0, and a bunch of dependencies. The install will fail with various error messages until you install those dependencies.
library(MIMOSA) #load MIMOSA
```

### Running the model

```r
library(MIMOSA)
data<-read.csv("mydata.csv")
E<-ConstructMIMOSAExpressionSet(data,
  reference=STIMULATION%in%"Unstimulated",
  measure.columns=c("NSUB","CYTNUM"),
  .variables=.(SUBJECTID,VISIT,CYTOKINE,TCELL)) #There's a lot of options here. See the documentation and vignette.
MIMOSA(NSUB+CYTNUM~SUBJECTID|TCELL+VISIT,
  E,
  subset=RefTreat%in%"Treatment"&CYTOKINE%in%"IL2",
  ref=RefTreat%in%"Reference"&CYTOKINE%in%"IL2")
```

A few things are worth explaining in the code above.  
ConstructMIMOSAExpressionSet takes a data.frame as input.
The expected information in the data frame is described in the documentation and the vignette.  
- It should contain positive and negative cell counts.  
  - Above, it's assumed they're in the *NSUB* and *CYTNUM* columns.  
- There should be a column describing the stimulation applied to each sample.  
  - Above, it's assumed that it's in the *STIMULATION* column. 
- The *reference* argument tells the code which stimulation is the reference. In this case it's the samples labelled 'Unstimulated'. 
- The *.variables* argument describes the grouping variables used to stratify the data into a single experimental unit.
  - In this imaginary data set, we have stimulated and unstimualted (*STIMULATION*) samples from different subjects (*PTID*), and at different timepoints (*VISIT*), measuring different responses (*CYTOKINE*), and different cells (*TCELL*).
  - We stratify by *PTID*, *VISIT*, *CYTOKINE*, and *TCELL*. Within each unique set of levels in that group, we'll have all the different stimulations for a subject. 
  - The code reshapes the data so that each stimulated sample is paired with the appropriate unstimulated control (defined in the *reference* argument). That is why we don't include the *STIMULATION* in the *.variables* argument.
- The *other.annotations* argument (not shown) lets you specify which columns to include in the phenodata for the expression set.

The call to *MIMOSA* takes the ExpressionSet as input. It supports the formula interface.  
The formula should have the negative and positive cell counts on the left hand side (in this case NSUB is the negative cell count, and CYTNUM is the positive cell count). The right hand describes how to fit the model. In this case we fit a separate model for each visit and T cell subset across all subjects. We restrict ourselves to the IL2 cytokine. RefTreat is a pre-defined variable with levels "Reference" and "Treatment" defining the stimulated and unstimualted groups.  
The definition of the reference and treatment groups will be handled internally and hidden from the end user in a future update.
