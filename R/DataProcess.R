
#############################################
## dataProcess
#############################################
#' @export dataProcess
#' @import survival 
#' @import preprocessCore 
#' @import MASS
#' @importFrom reshape2 dcast melt
#' @importFrom stats medpolish aggregate t.test lm summary.lm fitted resid p.adjust
#' @importFrom stats C approx coef cor dist formula loess median na.omit
#' @importFrom stats predict pt qnorm qt quantile reshape rnorm runif sd var vcov xtabs
#' @importFrom utils head read.table sessionInfo write.csv write.table
#' @importFrom methods validObject
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr filter n
#' @importFrom progress progress_bar
#' @importFrom tidyr gather

dataProcess  <-  function(raw,
                          logTrans=2,
                          normalization="equalizeMedians",
                          nameStandards=NULL,
                          address="",
                          fillIncompleteRows=TRUE,
                          featureSubset="all",
                          remove_noninformative_feature_outlier=FALSE,
                          n_top_feature=3,
                          summaryMethod="TMP",
                          equalFeatureVar=TRUE,
                          censoredInt="NA",
                          cutoffCensored="minFeature",
                          MBimpute=TRUE,
                          remove50missing=FALSE,
                          maxQuantileforCensored=0.999) {
  
    ## save process output in each step
    allfiles <- list.files()
  
    num <- 0
    filenaming <- "msstats"
    finalfile <- "msstats.log"
  
    while(is.element(finalfile, allfiles)) {
        num <- num + 1
        finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
    }
  
    session <- sessionInfo()
    sink("sessionInfo.txt")
    print(session)
    sink()
  
    processout <- as.matrix(read.table("sessionInfo.txt", header=TRUE, sep="\t"))
    write.table(processout, file=finalfile, row.names=FALSE)
  
    processout <- rbind(processout, as.matrix(c(" "," ","MSstats - dataProcess function"," "), ncol=1))
  
    ## make case-insensitive for function options
    ## ------------------------------------------
  
    normalization <- toupper(normalization)

    ## Check correct option or input
    ## check right column in input
  
    requiredinput <- c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                     "FragmentIon", "ProductCharge", 
                     "Condition", "BioReplicate", "Run", "Intensity")

    ## [THT: disambiguation for PeptideSequence & PeptideModifiedSequence - begin]
    ## PeptideModifiedSequence is also allowed.
    requiredInputUpper <- toupper(requiredinput)
    providedInputUpper <- toupper(colnames(raw))
    if (all(requiredInputUpper %in% providedInputUpper)) {
        processout <- rbind(processout, c("The required input : provided - okay"))
        write.table(processout, file = finalfile, row.names = FALSE)
    } else if (all(setdiff(requiredInputUpper, "PEPTIDESEQUENCE") %in% providedInputUpper) && "PEPTIDEMODIFIEDSEQUENCE" %in% providedInputUpper) {
        processout <- rbind(processout, c("The required input : provided - okay"))
        write.table(processout, file = finalfile, row.names = FALSE)
        # if PeptideModifiedSequence is provided instead of PeptideSequence, 
        # change the column name as PeptideSequence
        colnames(raw)[which(providedInputUpper == "PEPTIDEMODIFIEDSEQUENCE")]  <-  "PeptideSequence"
    } else {
        missedInput <- which(!(requiredInputUpper %in% providedInputUpper))
        processout <- rbind(processout, c(paste("ERROR : The required input : ", 
                            paste(requiredinput[missedInput], collapse = ", "), 
                            " are not provided in input - stop")))
        write.table(processout, file = finalfile, row.names = FALSE)
        stop("Please check the required input. The required input needs (ProteinName, PeptideSequence (or PeptideModifiedSequence), PrecursorCharge, FragmentIon, ProductCharge, Condition, BioReplicate, Run, Intensity)")
    }
    ## [THT: disambiguation for PeptideSequence & PeptideModifiedSequence - end]
  
    
    ## check logTrans is 2,10 or not
    if (logTrans!=2 & logTrans!=10) {
        processout <- rbind(processout,c("ERROR : Logarithm transformation : log2 or log10 only - stop"))
        write.table(processout, file=finalfile,row.names=FALSE)
      
        stop("Only log2 or log10 are possible.\n")
    }
    
    ## check no row for some feature : balanced structure or not  
    if (!(fillIncompleteRows==TRUE | fillIncompleteRows==FALSE) | !is.logical(fillIncompleteRows)) {
        processout <- rbind(processout, c(paste("The required input - fillIncompleteRows : 'fillIncompleteRows' value is wrong. It should be either TRUE or FALSE. - stop")))
        write.table(processout, file=finalfile, row.names=FALSE)
      
      stop("'fillIncompleteRows' must be one of TRUE or FALSE as a logical value.")
    }
    
    ## check input for summaryMethod
    
    if (sum(summaryMethod == c("linear", "TMP")) == 0) {
        processout <- rbind(processout,c("The required input - summaryMethod : 'summaryMethod' value is wrong. It should be one of 'TMP' or 'linear'. - stop"))
        write.table(processout, file=finalfile, row.names=FALSE)
      
        stop("'summaryMethod' value is wrong. It should be one of 'TMP' or 'linear'.")
    } else {
        processout <- rbind(processout, c(paste("summaryMethod : ", as.character(summaryMethod), sep="")))
        write.table(processout, file=finalfile, row.names=FALSE)
    }
    
    ## check input for cutoffCensored
    if (sum(cutoffCensored==c("minFeature","minRun","minFeatureNRun"))==0) {
        processout <- rbind(processout,c("The required input - cutoffCensored : 'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'. - stop"))
        write.table(processout, file=finalfile, row.names=FALSE)
      
        stop("'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'.")
    } else {
        processout <- rbind(processout,c(paste("cutoffCensored : ",as.character(cutoffCensored), sep="")))
        write.table(processout, file=finalfile, row.names=FALSE)
    }
    
    ## check input for censoredInt
    if (sum(censoredInt == c("0", "NA")) == 0 & !is.null(censoredInt)) {
        processout <- rbind(processout,c("The required input - censoredInt : 'censoredInt' value is wrong. 
                                         It should be one of '0','NA', NULL. - stop"))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("'censoredInt' value is wrong. It should be one of '0','NA',NULL.")
    } else {
        processout <- rbind(processout, c(paste("censoredInt : ", as.character(censoredInt), sep="")))
        write.table(processout, file=finalfile, row.names=FALSE)
    }
    
    ## check input for censoredInt and MBimpute
    if ( summaryMethod == 'TMP' & MBimpute & is.null(censoredInt) ) {
        processout <- rbind(processout, c("The rcombination of equired input - 
                                          censoredInt and MBimpute : 'censoredInt=NULL' has no censored missing values. 
                                          Imputation will not be performed.- stop"))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("'censoredInt=NULL' means that dataset has no censored missing value and MSstats will not impute.
             But, 'MBimpute=TRUE' is selected. Please replace by 'MBimpute=FALSE' or censoredInt='NA' or '0'")
      
    } 
    
    ## [THT: if (!all(normalization %in% c("NONE", "FALSE", "EQUALIZEMEDIANS", "QUANTILE", "GLOBALSTANDARDS")))]
    ## [THT: send a warning message if the user mixes "NONE" with any of the last three choices]
    if (!(normalization=="NONE" | normalization=="FALSE" | 
          normalization=="EQUALIZEMEDIANS" | normalization=="QUANTILE" | 
          normalization=="GLOBALSTANDARDS")) {
        processout <- rbind(processout,c(paste("The required input - normalization : 'normalization' value is wrong. - stop")))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("'normalization' must be one of \"None\", \"FALSE\", \"equalizeMedians\", 
             \"quantile\", or \"globalStandards\". Please assign 'normalization' again.")
    } 
    
    ## need the names of global standards
    if (!is.element("NONE",normalization) & 
        !is.element("FALSE",normalization) & 
        is.element("GLOBALSTANDARDS",normalization) & 
        is.null(nameStandards)) {
        processout <- rbind(processout, c("ERROR : For normalization with global standards, 
                                          the names of global standards are needed. Please add 'nameStandards' input."))
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop ("For normalization with global standards, the names of global standards are needed. 
              Please add 'nameStandards' input." )
    }
    
    ## check whether class of intensity is factor or chaterer, if yes, neec to chage as numeric
    if (is.factor(raw$Intensity) | is.character(raw$Intensity)) {   
        suppressWarnings(raw$Intensity <- as.numeric(as.character(raw$Intensity)))
    }
  
    ## check whether the intensity has 0 value or negative value
    #   if (length(which(raw$Intensity<=0))>0 & !skylineReport) {
    
    #   if (is.null(censoredInt)) {
    #       processout <- rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
    #       write.table(processout, file=finalfile,row.names=FALSE)
    
    #       stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
        
    #   } else if (censoredInt=="NA") {
            
    #       processout <- rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
    #       write.table(processout, file=finalfile,row.names=FALSE)
    
    #       stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
            
  #     }
    #}
    

    ## here, need to get standard protein name
    ## column name : standardtype..
    ## what value it has, normzalition, unique(proteinname)
    ## if normalition== "standard" & no normalizaion selection, error message
  
    ## annotation information : 
    if ( any(is.na(raw$Run)) ) {
        processout <- rbind(processout, c("ERROR : There is missing information in 'Run' column. Please check 'Run' column."))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop ("There is missing information in 'Run' column. Please check 'Run' column." )
    }
    
    if ( any(is.na(raw$BioReplicate)) ) {
        processout <- rbind(processout, c("ERROR : There is missing information in 'BioReplicate' column. 
                                          Please check 'BioReplicate' column."))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop ("There is missing information in 'BioReplicate' column. Please check 'BioReplicate' column." )
    }
    
    if ( any(is.na(raw$Condition)) ) {
        processout <- rbind(processout, c("ERROR : There is missing information in 'Condition' column. 
                                          Please check 'Condition' column."))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop ("There is missing information in 'Condition' column. Please check 'Condition' column." )
    }
    
    ## make letters case-insensitive
    colnames(raw) <- toupper(colnames(raw))
    
    require.col <- c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", 
                     "FRAGMENTION", "PRODUCTCHARGE", 
                     "CONDITION", "BIOREPLICATE", "RUN", "INTENSITY")
    raw.temp <- raw[, require.col]
  
    ## before remove, get PeptideSequence and combination of PeptideSequence and precursorcharge for global standard normalization
    tempPeptide <- unique(raw[, c("PEPTIDESEQUENCE", "PRECURSORCHARGE")])
    tempPeptide$PEPTIDE <- paste(tempPeptide$PEPTIDESEQUENCE, tempPeptide$PRECURSORCHARGE, sep="_")
  
    rm(raw)
  
    ## assign peptide, transition
    raw.temp <- data.frame(raw.temp, 
                           PEPTIDE=paste(raw.temp$PEPTIDESEQUENCE, raw.temp$PRECURSORCHARGE, sep="_"), 
                           TRANSITION=paste(raw.temp$FRAGMENTION, raw.temp$PRODUCTCHARGE, sep="_"))
  
    require.col <- c("PROTEINNAME", "PEPTIDE", "TRANSITION", 
                     "CONDITION", "BIOREPLICATE", "RUN",  "INTENSITY")
    
    raw.temp <- raw.temp[, require.col]
  
    colnames(raw.temp) <- c("Protein", "Peptide", "Transition", 
                            "Condition", "Sample", "Run", "Intensity")
 
    ## create work data for quant analysis
    ## -----------------------------------
    raw.temp <- raw.temp[!is.na(raw.temp$Protein), ]
    raw.temp <- raw.temp[raw.temp$Protein != '', ]
    
    work <- data.frame(PROTEIN=raw.temp$Protein, 
                       PEPTIDE=raw.temp$Peptide, 
                       TRANSITION=raw.temp$Transition, 
                       FEATURE=paste(raw.temp$Peptide, raw.temp$Transition, sep="_"), 
                       GROUP_ORIGINAL=raw.temp$Condition, 
                       SUBJECT_ORIGINAL=raw.temp$Sample, 
                       RUN=raw.temp$Run, 
                       GROUP=0,
                       SUBJECT=0,
                       INTENSITY=raw.temp$Intensity)
  
    work$GROUP_ORIGINAL <- factor(work$GROUP_ORIGINAL)
    work$SUBJECT_ORIGINAL <- factor(work$SUBJECT_ORIGINAL, levels=unique(work$SUBJECT_ORIGINAL))

    work[, "GROUP"] <- work[, "GROUP_ORIGINAL"]
    work[, "SUBJECT"] <- work[, "SUBJECT_ORIGINAL"]

    work <- data.frame(work, SUBJECT_NESTED=paste(work$GROUP, work$SUBJECT, sep="."))
    
    processout <- rbind(processout, c("New input format : made new columns for analysis - okay"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
  
    ## 2016. 08.29 : replace <1 with zero for log2(intensity)
    if ( length(which(!is.na(work$INTENSITY) & work$INTENSITY < 1)) > 0 ) {
      
      processout <- rbind(processout, c(paste0("** There are ",  
                                              length(which(!is.na(work$INTENSITY) & work$INTENSITY < 1)), 
                                              " intensities which are zero. These intensities are replaced with 1.")))
      write.table(processout, file=finalfile, row.names=FALSE)
      
      message(paste0("** There are ", length(which(!is.na(work$INTENSITY) & work$INTENSITY < 1)), 
                    " intensities which are zero or less than 1. These intensities are replaced with 1."))
      
      work[!is.na(work$INTENSITY) & work$INTENSITY < 1, 'INTENSITY'] <- 1
    }
    
    ## log transformation
    work$ABUNDANCE <- work$INTENSITY
  
    ## now, INTENSITY keeps original values.
    
    ## NA means no observation. assume that spectral tools are not report if no observation. zero means detected but zero. 
    ## considered intenseity <1 -> intensity = 1
    ## work[!is.na(work$ABUNDANCE) & work$ABUNDANCE==0,"ABUNDANCE"] <- 1
    
    ## based on logTrans option, assign log transformation
    ## remove log2 or log10 intensity
    ### [THT: add one more conidtion to have the program complain if a user 
    ### provide unexpected value for logTrans]
    if (logTrans == 2) {
        work$ABUNDANCE <- log2(work$ABUNDANCE)
    } else if (logTrans == 10) {
        work$ABUNDANCE <- log10(work$ABUNDANCE)
    }   
  
    processout <- rbind(processout, 
                        c(paste0("Logarithm transformation: log", logTrans, 
                                " transformation is done - okay")))
    write.table(processout, file=finalfile, row.names=FALSE)
   
    work$RUN <- factor(work$RUN)

    ## check messingness for multirun 
  
    ## check no value for some feature : balanced structure or not
    ## need to separate label-free or label-based
  
    processout <- rbind(processout, c(paste("fillIncompleteRows = ", fillIncompleteRows,sep="")))
    write.table(processout, file=finalfile, row.names=FALSE)
  
    ## [THT: better to write a function for single method, and call that function
    ## here and for the case with multuple methods]
    ## only 1 method
  
    ## label-free experiments

    ## get feature by Run count of data
    structure = tapply ( work$ABUNDANCE, list ( work$FEATURE, work$RUN ) , function ( x ) length ( x ) ) 

    ## structure value should be 1 for label-free, if not there are missingness. if more there are duplicates.

    flagmissing = sum(is.na(structure)) > 0
    flagduplicate = sum(structure[!is.na(structure)]>1) > 0

    ### if there is missing rows
    if ( flagmissing ) {
        processout <- rbind(processout, c("CAUTION: the input dataset has incomplete rows. 
                                          If missing peaks occur they should be included in the dataset as separate rows, 
                                          and the missing intensity values should be indicated with 'NA'. 
                                          The incomplete rows are listed below."))
        write.table(processout, file=finalfile,row.names=FALSE)

        message("CAUTION : the input dataset has incomplete rows. 
                If missing peaks occur they should be included in the dataset as separate rows, 
                and the missing intensity values should be indicated with 'NA'. 
                The incomplete rows are listed below.")

        ## first, which run has missing 
        runstructure <- apply ( structure, 2, function ( x ) sum ( is.na ( x ) ) ) > 0

        ## get the name of Run
        runID <- names(runstructure[runstructure==TRUE])

        ## for missign row, need to assign before looping
        missingwork <- NULL

        ## then for each run, which features are missing,
        for(j in 1:length(runID)) {
  
            ## get subject, group information for this run
            nameID <- unique(work[work$RUN==runID[j], c("SUBJECT_ORIGINAL","GROUP_ORIGINAL",
                                                            "GROUP","SUBJECT","SUBJECT_NESTED",
                                                            "RUN")])
  
            ## get feature ID
            featureID <- structure[,colnames(structure)==runID[j]]
  
            ## get feature ID which has no measuremnt.
            finalfeatureID <- featureID[is.na(featureID)]
  
            ## print features ID        
            message(paste0("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]),
                              ", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), 
                              " has incomplete rows for some features (", 
                              paste(names(finalfeatureID), collapse=", "), ")"))
  
            ## save in process file.
            processout <- rbind(processout, c(paste0("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]),
                                                        ", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), 
                                                        " has incomplete rows for some features (", 
                                                        paste(names(featureID[is.na(featureID)]), collapse=", "), ")")))
            write.table(processout, file=finalfile, row.names=FALSE)
  
            ## add missing rows if option is TRUE
            if (fillIncompleteRows) {
    
                tempTogetfeature <- work[which(work$FEATURE %in% names(finalfeatureID)), ]
    
                ## get PROTEIN and FEATURE infomation
                tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
    
                ## merge feature info and run info as 'work' format
                tempmissingwork <- data.frame(tempfeatureID, 
                                              GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, 
                                              SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, 
                                              RUN=nameID$RUN, 
                                              GROUP=nameID$GROUP, 
                                              SUBJECT=nameID$SUBJECT, 
                                              SUBJECT_NESTED=nameID$SUBJECT_NESTED, 
                                              INTENSITY=NA, 
                                              ABUNDANCE=NA) 
    
                ## merge with tempary space, missingwork
                missingwork <- rbind(missingwork, tempmissingwork)
            } # end fillIncompleteRows options
        } # end loop for run ID

        ## [THT: this part can probably be merged into the above. 
        ## Also, it might be better to check fillIncompleteRows earlier
        ## and terminate the process when it's FALSE]
        if (fillIncompleteRows) {
  
            ## merge with work
            ## in future, use rbindlist?? rbindlist(list(work, missingwork))
            work <- rbind(work, missingwork)
  
            ## print message
            message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
  
            ## save in process file.
            processout <- rbind(processout, "Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
            write.table(processout, file=finalfile, row.names=FALSE)
  
        } else {
  
            ## save in process file.
            processout <- rbind(processout,"Please check whether features in the list are generated from spectral processing tool. 
                                Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
            write.table(processout, file=finalfile,row.names=FALSE)
  
            stop("Please check whether features in the list are generated from spectral processing tool or not. 
                 Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
  
        }
    } # end for flag missing
    
    ## if there are duplicates measurements
    if (flagduplicate) {

        ## first, which run has duplicates
        runstructure <- apply ( structure, 2, function ( x ) sum (x[!is.na(x)] > 1 ) > 0 )

        runID <- names(runstructure[runstructure==TRUE])

        ## then for each run, which features have duplicates,
        for(j in 1:length(runID)) {
  
            nameID <- unique(work[work$RUN == runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", 
                                                          "GROUP","SUBJECT", "SUBJECT_NESTED", 
                                                          "RUN")])
  
            featureID <- structure[, colnames(structure)==runID[j]]
            finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
  
            message(paste0("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]),
                          ", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), 
                          " has multiple rows (duplicate rows) for some features (", 
                          paste(names(finalfeatureID), collapse=", "), ")"))
  
            ## save in process file.
            processout <- rbind(processout, c(paste0("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]), 
                                                    ", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), 
                                                    " has multiple rows (duplicate rows) for some features (", 
                                                    paste(names(featureID[is.na(featureID)]), collapse=", "), ")")))
            write.table(processout, file=finalfile, row.names=FALSE)
        }

        ## save in process file.
        processout <- rbind(processout,"Please remove duplicate rows in the list above. ")
        write.table(processout, file=finalfile,row.names=FALSE)

        stop("Please remove duplicate rows in the list above.\n")       
    } # end flag duplicate

        ## no missing and no duplicates
    if (!flagmissing & !flagduplicate) {
        processout <- rbind(processout, c("Balanced data format with NA for missing feature intensities - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
    } 
    ## end label-free

    ## factorize GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN
    ## -------------------------------------------------------------------------------------------------

    work$PROTEIN <- factor(work$PROTEIN)
    work$PEPTIDE <- factor(work$PEPTIDE)
    work$TRANSITION <- factor(work$TRANSITION)
  
    work <- work[with(work, order(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, PROTEIN, PEPTIDE, TRANSITION)),]
  
    work$GROUP <- factor(work$GROUP)
    work$SUBJECT <- factor(work$SUBJECT)
    ## SUBJECT_ORIGINAL_NESTED will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL
  
    work$SUBJECT_NESTED <- factor(work$SUBJECT_NESTED, levels=unique(work$SUBJECT_NESTED))
  
    ## FEATURE will sorted as PROTEIN, PEPTIDE, TRANSITION
    work$FEATURE <- factor(work$FEATURE, levels=unique(work$FEATURE))

    ## RUN will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN
    work$originalRUN <- work$RUN
    work$RUN <- factor(work$RUN, levels=unique(work$RUN), labels=seq(1, length(unique(work$RUN))))
  
    processout <- rbind(processout, c("Factorize in columns(GROUP, SUBJECT, GROUP_ORIGINAL, 
                                      SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN) - okay"))
    write.table(processout, file=finalfile, row.names=FALSE)
  
  
    ## Normalization ##
    ## ------------- ##
       
    ## Normalization : option 0. none
    if (is.element("NONE",normalization) | is.element("FALSE",normalization)) { # after 'toupper', FALSE becomes character.
        processout <- rbind(processout, c("Normalization : no normalization - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
    }
  

    ## Normalization : option 1. constant normalization , equalize medians ##
    ## -------------------------------------------------------------------
    if (!is.element("NONE", normalization) & 
        !is.element("FALSE", normalization) & 
        is.element("EQUALIZEMEDIANS", normalization)) {
    
        ## Constant normalization by endogenous per method
  
        ## [MC : use median of medians]
        median.run.method  <-  aggregate(ABUNDANCE ~ RUN , data = work, median, na.rm = TRUE)
        global.median <- median(median.run.method$ABUNDANCE, na.rm = TRUE)
  
        namerun <- unique(work[, "RUN"])

        for (i in 1:length(namerun)) {
            ## ABUNDANCE is normalized
          namerun.idx <- which(work$RUN == namerun[i])
          work[namerun.idx, "ABUNDANCE"] <- work[namerun.idx, "ABUNDANCE"] - median.run.method[median.run.method$RUN == namerun[i], "ABUNDANCE"] + global.median
        }
    
        processout <- rbind(processout, c("Normalization : Constant normalization (equalize medians) - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
    } ## end equaliemedian normalization    
  
    ## Normalization : option 2. quantile normalization ##
    ## ------------------------------------------------ ##
    if (!is.element("NONE", normalization) & 
        !is.element("FALSE", normalization) & 
        is.element("QUANTILE", normalization)) {
    
        ## for label-free, just use endogenous
  
        quantileall <- NULL
        
        ## ABUNDANCE=0 replace with 1, in order to distinguish later.
        work[!is.na(work$ABUNDANCE) & work$ABUNDANCE == 0, 'ABUNDANCE'] <- 1
  
        worktemp <- work[which(!is.na(work$INTENSITY)),]
        worktemp$RUN <- factor(worktemp$RUN)
        worktemp$FEATURE <- factor(worktemp$FEATURE)

        quantiletemp <- as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))

        ## need to put NA for missing value in endogenous
        quantiletemp[quantiletemp == 0] <- NA

        ## using preprocessCore library
        quantiledone <- normalize.quantiles(quantiletemp)
        rownames(quantiledone) <- rownames(quantiletemp)
        colnames(quantiledone) <- colnames(quantiletemp)
        
        ## get quantiled to long format for apply difference endogenous
        quantilelong <- melt(quantiledone, id=rownames(quantiledone))
        colnames(quantilelong) <- c("FEATURE", "RUN", "ABUNDANCE_quantile")
        rm(quantiledone)

        ## quantileall <- rbindlist(list(quantileall,quantilelong))
        quantileall <- rbind(quantileall, quantilelong)

        rm(quantilelong)
  
        work <- merge(work, quantileall, by=c("FEATURE", "RUN"))
        rm(quantileall)
  
        ## reorder
        work <- data.frame("PROTEIN"=work$PROTEIN, 
                           "PEPTIDE"=work$PEPTIDE, 
                           "TRANSITION"=work$TRANSITION, 
                           "FEATURE"=work$FEATURE, 
                           "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, 
                           "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, 
                           "RUN"=work$RUN, 
                           "GROUP"=work$GROUP, 
                           "SUBJECT"=work$SUBJECT, 
                           "SUBJECT_NESTED"=work$SUBJECT_NESTED, 
                           "INTENSITY"=work$INTENSITY, 
                           "ABUNDANCE"=work$ABUNDANCE_quantile, 
                           "originalRUN"=work$originalRUN)
  
        work <- work[with(work, order(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, PROTEIN, PEPTIDE, TRANSITION)), ]
      
        ## for skyline case, separate 1 and zero
        work[!is.na(work$INTENSITY) & work$INTENSITY == 1, 'ABUNDANCE'] <- 0
    
        processout <- rbind(processout, c("Normalization : Quantile normalization - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
    }
  
  
    ## Normalization : option 3. global standards - for endogenous ##
    ## ----------------------------------------------------------- ##
    if (!is.element("NONE", normalization) & 
        !is.element("FALSE", normalization) & 
        is.element("GLOBALSTANDARDS", normalization)) {
    
        work$RUN <- factor(work$RUN)
        combine <- data.frame(RUN=levels(work$RUN))
        allPeptide <- unique(work$PEPTIDE)
        allProtein <- unique(work$PROTEIN)
    
        for (i in 1:length(nameStandards)) {
      
            ## if Peptides
            ## namePeptide <- allPeptide[grep(nameStandards[i],allPeptide)] ## cannot grep for modified peptide sequence, [,],+ sign
            namePeptide <- tempPeptide[tempPeptide$PEPTIDESEQUENCE == nameStandards[i], "PEPTIDE"]
      
            if (length(namePeptide)!=0) {
                tempStandard <- work[work$PEPTIDE == namePeptide,]
            } else {
        
                ## if Proteins
                nameProtein <- allProtein[allProtein == nameStandards[i]] # if we use 'grep', can' find the proteins name with some symbol, such as 'sp|P30153|2AAA_HUMAN'
        
                if (length(nameProtein)!=0) {
                    tempStandard <- work[work$PROTEIN==nameProtein,]
                } else {
                    processout <- rbind(processout,c(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not.")))
                    write.table(processout, file=finalfile,row.names=FALSE)
          
                    stop(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not."))
                }   
            }
      
            ## here, by RUN, but need to check !!!
            tempStandard <- tempStandard[tempStandard$GROUP!="0",]
            tempStandard$RUN <- factor(tempStandard$RUN)
      
            tempStandard <- tempStandard[!is.na(tempStandard$ABUNDANCE),]
            meanStandard <- tapply(tempStandard$ABUNDANCE, tempStandard$RUN, function(x) mean(x, na.rm=TRUE))
      
            meanStandard <- data.frame(RUN=names(meanStandard),meanStandard)
            combine <- merge(combine, meanStandard, by="RUN", all=TRUE)
            colnames(combine)[i+1] <- paste("meanStandard",i,sep="")
        }
    
        rownames(combine) <- combine$RUN
        combine <- subset(combine, select=-c(RUN))
    
        ## get mean among global standards
        allmean <- apply(combine,1, function(x) mean(x, na.rm=TRUE))
        ## allmean[is.na(allmean)] <- 0
	median.all <- median(allmeantemp$allmean, na.rm=TRUE)
    
        ## adjust
        namerun <- unique(work[, "RUN"])
  
        for (i in 1:length(namerun)) {
            ## ABUNDANCE is normalized          
            if (!is.na(allmean[names(allmean)==namerun[i]])) work[work$RUN==namerun[i],"ABUNDANCE"] <- work[work$RUN==namerun[i],"ABUNDANCE"]-allmean[names(allmean)==namerun[i]] + median.all
        }
    
        processout <- rbind(processout, c("Normalization : normalization with global standards protein - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
    
    }
    
    ## ----------------------------------------------------------- ##
    ## after normalization, zero intensity could be negative
    
    ## if abundance became less than zero, after normalization
    work[!is.na(work$ABUNDANCE) & work$ABUNDANCE < 0, "ABUNDANCE"] <- 0
    
    ## if abundance become greater than zero, after normalization. 
    ## hard to know how much higher, so, use intensity value, which is not used for noramlization
    work[!is.na(work$INTENSITY) & work$INTENSITY == 1, "ABUNDANCE"] <- 0
    
    #Below two lines were merely for in-house testing and comparisons when needed
    #work.NoImpute <- work
    #AbundanceAfterImpute <- .Imputation(work, cutoffCensored, censoredInt, remove50missing, MBimpute, original_scale)
    
    ## ------------- ## 
    ## how to decide censored or not
    ## ------------- ##
    
    ### If imputation=TRUE and there is any value for maxQuantileforCensored, apply cutoff for censored missing
    if ( summaryMethod == "TMP" & MBimpute ) {
      
        work$censored <- FALSE
        
        ## if intensity = 1, but abundance > cutoff after normalization, it also should be censored.
        
        if( !is.null(maxQuantileforCensored) ) {
            ### calculate outlier cutoff
            ## only consider intensity > 1
            tmp <- work[!is.na(work$INTENSITY) & work$INTENSITY > 1, 'ABUNDANCE']
            ## or
            #tmp <- work[!is.na(work$INTENSITY), 'ABUNDANCE']
    
            log2int.prime.quant <- quantile(tmp, prob=c(0.01, 0.25, 0.5, 0.75, maxQuantileforCensored), na.rm = TRUE)
            iqr <- log2int.prime.quant[4] - log2int.prime.quant[2]
    
            ### need to decide the multiplier from high intensities
            multiplier <- (log2int.prime.quant[5] - log2int.prime.quant[4])/iqr
    
            cutoff.lower <- (log2int.prime.quant[2] - multiplier * iqr) 
    
            work[!is.na(work$INTENSITY) & 
                     work$ABUNDANCE < cutoff.lower, 'censored'] <- TRUE
            
            message(paste('** Log2 intensities under cutoff =', 
                          format(cutoff.lower, digits=5), 
                          ' were considered as censored missing values.'))
            
            processout <- rbind(processout, 
                                c(paste('** Log2 intensities under cutoff =', 
                                        format(cutoff.lower, digits=5), 
                                        ' were considered as censored missing values.')))
            write.table(processout, file=finalfile, row.names=FALSE)
            
            ## if censoredInt == '0, and cutoff is negative, still zero should becensored
            if ( cutoff.lower <= 0 & !is.null(censoredInt) & censoredInt == "0" ) {
                
                work[!is.na(work$INTENSITY) & work$INTENSITY == 1, 'censored'] <- TRUE
                work[!is.na(work$ABUNDANCE) & work$ABUNDANCE <= 0, 'censored'] <- TRUE
                
                message(paste('** Log2 intensities = 0 were considered as censored missing values.'))
                
                processout <- rbind(processout, 
                                    c(paste('** Log2 intensities = 0 were considered as censored missing values.')))
                write.table(processout, file=finalfile, row.names=FALSE)
                
            }
            
            ## if censoredInt == NA, original NA also shoule be 'censored'
            if (!is.null(censoredInt) & censoredInt == "NA") {
      
                work[is.na(work$INTENSITY), 'censored'] <- TRUE
                
                message(paste('** Log2 intensities = NA were considered as censored missing values.'))
                
                processout <- rbind(processout, c('** Log2 intensities = NA were considered as censored missing values.'))
                write.table(processout, file=finalfile, row.names=FALSE)
      
            }
      
        } else { ## will MBimpute, but not apply algorithm for cutoff
        
            if(censoredInt == '0'){
                work[!is.na(work$INTENSITY) & work$INTENSITY == 1, 'censored'] <- TRUE
                work[ !is.na(work$ABUNDANCE) & work$ABUNDANCE <= 0, 'censored'] <- TRUE
            }
            if(censoredInt == 'NA'){
                work[is.na(work$ABUNDANCE), 'censored'] <- TRUE
            }
        
        }
      
    }
    

    ## ------------- ## 
    ## featureSubset ##
    ## ------------- ##
    ##  !! need to decide how to present : keep original all data and make new column to mark, or just present selected subset    
  
    if (featureSubset == "all") {
        message("** Use all features that the dataset origianally has.")
     
        processout <- rbind(processout, c("** Use all features that the dataset origianally has."))
        write.table(processout, file=finalfile, row.names=FALSE)
    } 

    if (featureSubset == "highQuality") {
        
        message("** Selecting high quality features temporarily defaults to featureSubset = top3. 
                Updates for this option will be available in the next release.")
        
        featureSubset <- 'top3'
        
        processout <- rbind(processout, c("** Selecting high quality features temporarily defaults to featureSubset = top3. 
                                          Updates for this option will be available in the next release."))
        
        write.table(processout, file=finalfile, row.names=FALSE)
    }
  
    if (featureSubset == "top3") {
        message("** Use top3 features that have highest average of log2(intensity) across runs.")
        
        processout <- rbind(processout, c("** Use top3 features that have highest average of log2(intensity) across runs."))
        write.table(processout, file=finalfile, row.names=FALSE)
     
        ## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
        ## how to decide top3 for DIA?
        work$remove <- FALSE
    
        temp1 <- aggregate(INTENSITY~PROTEIN+FEATURE,data=work, function(x) mean(x, na.rm=TRUE))

        temp2 <- split(temp1, temp1$PROTEIN)

        temp3 <- lapply(temp2, function(x) { 
            x <- x[order(x$INTENSITY, decreasing=TRUE), ]
            x <- x$FEATURE[1:3]
            })
    
        selectfeature <- unlist(temp3, use.names=FALSE)
        selectfeature <- selectfeature[!is.na(selectfeature)]
        
        ## get subset
        work[-which(work$FEATURE %in% selectfeature), 'remove'] <- TRUE

    }
  
    if (featureSubset == "topN") {
    
        ## check whether there is the input for 'N'
    
        message(paste0("** Use top", n_top_feature, " features that have highest average of log2(intensity) across runs."))
      
        processout <- rbind(processout, c(paste0("** Use top", n_top_feature, 
                                                " features that have highest average of log2(intensity) across runs.")))
        write.table(processout, file=finalfile, row.names=FALSE)
      
        ## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
        ## how to decide top3 for DIA?
      
        work$remove <- FALSE
      
        worktemp <- work[!is.na(work$ABUNDANCE) & work$ABUNDANCE != 0, ]
        temp1 <- aggregate(INTENSITY ~ PROTEIN+FEATURE, data=worktemp, function(x) mean(x, na.rm=TRUE))
      
        temp2 <- split(temp1, temp1$PROTEIN)
      
        temp3 <- lapply(temp2, function(x) { 
            x <- x[order(x$INTENSITY, decreasing=TRUE), ]
            x <- x$FEATURE[1:n_top_feature]
        })
      
        selectfeature <- unlist(temp3, use.names=FALSE)
        selectfeature <- selectfeature[!is.na(selectfeature)]
      
        ## get subset
        work[-which(work$FEATURE %in% selectfeature), 'remove'] <- TRUE
      
    }
  
    ## check missingness 
    ## transitions are completely missing in one condition : missingness ##
    #Use the data frame before imputation to summarize the missingness
    all.work <- work    
    test <- tapply(is.na(work[, "ABUNDANCE"]), work[, c("GROUP_ORIGINAL", "FEATURE")], function(x) sum(x, na.rm=TRUE))
    numObs <- tapply(work[, "ABUNDANCE"], work[, c("GROUP_ORIGINAL", "FEATURE")], function(x) length(x))
    test1 <- test == numObs
    test2 <- apply(test1, 2, function(x) sum(x, na.rm=TRUE))
    filterList <- names(test2)[test2 > 0]
    final.decision <- ifelse(test2>0, 1, 0)
  
    ## output : summary ##
    ## ---------------- ##

    temp <- data.frame("Summary of Features :")
    colnames(temp) <- " "
    rownames(temp) <- " "
    print(temp)
  
    summary.f <- matrix(NA,nrow=3)
    summary.f[1] <- nlevels(work$PROTEIN)
  
    temp <- unique(work[, c("PROTEIN", "PEPTIDE")])
    temp1 <- xtabs(~PROTEIN, data=temp)
    temp2 <- summary(as.numeric(temp1))
    summary.f[2] <- paste(temp2["Min."], temp2["Max."], sep="-")
  
    temp <- unique(work[, c("PEPTIDE", "FEATURE")])
    temp1 <- xtabs(~PEPTIDE, data=temp)
    temp2 <- summary(as.numeric(temp1))
    summary.f[3] <- paste(temp2["Min."], temp2["Max."], sep="-")
  
    colnames(summary.f) <- "count"
    rownames(summary.f) <- c("# of Protein", "# of Peptides/Protein", "# of Transitions/Peptide")
  
    print(as.data.frame(summary.f))
  
    ## output for process
    processout <- rbind(processout, c("Summary of Features :"))
    processout <- rbind(processout, c(paste(rownames(summary.f)[1]," : ", summary.f[1], sep="")))
    processout <- rbind(processout, c(paste(rownames(summary.f)[2]," : ", summary.f[2], sep="")))
    processout <- rbind(processout, c(paste(rownames(summary.f)[3]," : ", summary.f[3], sep="")))
  
    write.table(processout, file=finalfile, row.names=FALSE)
  
    ## protein list with 1 feature
    temp <- unique(work[, c("PROTEIN", "FEATURE")])
    temp1 <- xtabs(~PROTEIN, data=temp)
    temp2 <- as.data.frame(temp1[temp1 == 1])

    if (nrow(temp2) > 0) {
        if(nrow(temp2) > 1){
            npro <- min(c(nrow(temp2), 10)) 
            message("\n","** " , nrow(temp2), 
                    " Proteins have only single transition : Consider excluding this protein from the dataset. (", 
                    paste(temp2$PROTEIN[1:npro], collapse = ", "), " ...) \n")  
        } else {
            message("\n","** " , nrow(temp2), 
                    " Proteins have only single transition : Consider excluding this protein from the dataset. (", 
                    rownames(temp2), ") \n")  
        }
       
    }
  
    temp <- data.frame("Summary of Samples :")
    colnames(temp) <- " "
    rownames(temp) <- " "
    print(temp)
  
    summary.s <- matrix(NA,ncol=nlevels(work$GROUP_ORIGINAL), nrow=3)
  
    ## # of MS runs
    temp <- unique(work[, c("GROUP_ORIGINAL", "RUN")])
    temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
    summary.s[1,] <- temp1
  
    ## # of biological replicates
    temp <- unique(work[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
    temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
    summary.s[2,] <- temp1
  
    ## # of technical replicates
    c.tech <- round(summary.s[1,] / summary.s[2,])
    ##summary.s[3,] <- ifelse(c.tech==1,0,c.tech)
    summary.s[3,] <- c.tech
  
    colnames(summary.s) <- unique(work$GROUP_ORIGINAL)
    rownames(summary.s) <- c("# of MS runs","# of Biological Replicates", "# of Technical Replicates")
  
    print(summary.s)
  
    message("\n Summary of Missingness :\n" )
    message("  # transitions are completely missing in one condition: ", sum(final.decision!=0), "\n")
    if (sum(final.decision!=0)!=0) {
        tmp.final <- final.decision[final.decision != 0]
        if( length(tmp.final) > 5 ){
            message("    -> ", paste(names(tmp.final[1:5]),collapse = ", "), " ...")
        } else {
            message("    -> ", paste(names(tmp.final),collapse = ", "), " ...")
        }
        rm(tmp.final)
    }
  
    without <- xtabs(~RUN, work)
    withall <- xtabs(~RUN, all.work)
    run.missing <- without / withall
    message("\n  # run with 75% missing observations: ", sum(run.missing<0.25), "\n")
    if (sum(run.missing<0.25)!=0) {
        message("    -> ", paste("RUN", names(without[run.missing<0.25]), sep=" "))
    }
  
    ## output process
    processout <- rbind(processout, c("Summary of Missingness :"))
    processout <- rbind(processout, c(paste0("  # transitions are completely missing in one condition: ", 
                                            sum(final.decision!=0))))
    if (sum(final.decision!=0)!=0){
        
        tmp.final <- final.decision[final.decision != 0]
        if( length(tmp.final) > 5 ){
            processout <- rbind(processout,"    -> ", paste(names(tmp.final[1:5]), collapse = ", "), " ...")
            
        } else {
            processout <- rbind(processout,"    -> ", paste(names(tmp.final), collapse = ", "), " ...")
        }
        rm(tmp.final)
        
    } 
  
    processout <- rbind(processout, c(paste0("  # run with 75% missing observations: ", sum(run.missing < 0.25))))
    if (sum(run.missing<0.25)!=0) {
        processout <- rbind(processout, "    -> ", paste("RUN", names(without[run.missing < 0.25]), sep=" "))
    }
    write.table(processout, file=finalfile, row.names=FALSE)
  
    processout <- rbind(processout, c("Processing data for analysis is done. - okay"))
    write.table(processout, file=finalfile, row.names=FALSE)
  

    ## get the summarization per subplot (per RUN)   
    ## -------------------------------------------
    
    message("\n == Start the summarization per subplot...")

    rqresult <- try(.runQuantification(work, summaryMethod, equalFeatureVar, 
                                       cutoffCensored, censoredInt, remove50missing, MBimpute, 
                                       original_scale=FALSE, logsum=FALSE, featureSubset,
                                       message.show=FALSE), silent=TRUE)

    if (class(rqresult) == "try-error") {
        message("*** error : can't summarize per subplot with ", summaryMethod, ".")
     
        processout <- rbind(processout, c(paste0("error : can't summarize per subplot with ", summaryMethod, ".")))
        write.table(processout, file=finalfile, row.names=FALSE)
    
        rqall <- NULL
        rqmodelqc <- NULL
        workpred <- NULL
      
    } else {
      
        if (sum(is.element(colnames(rqresult$rqdata), "RUN")) == 0) {
            ## logsum is summarization per subject
            lab <- unique(work[, c("GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
            rqall <- merge(rqresult$rqdata, lab, by="SUBJECT_ORIGINAL")
            
        } else {
            lab <- unique(work[, c("RUN", "originalRUN", "GROUP", "GROUP_ORIGINAL", 
                                   "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
            rqall <- merge(rqresult$rqdata, lab, by="RUN")
        }
        
        rqall$GROUP <- factor(rqall$GROUP)
        rqall$Protein <- factor(rqall$Protein)
        
        rqmodelqc <- rqresult$ModelQC
        
        #MC : can't use this predicted value.
        #workpred <- rqresult$PredictedBySurvival
        workpred <- NULL
        
        message("\n == the summarization per subplot is done.")
        
        processout <- rbind(processout, c(paste0("the summarization per subplot is done.- okay : ", summaryMethod)))
        write.table(processout, file=finalfile, row.names=FALSE)

    }
     
    ## return work data.frame   and run quantification
    
    #Align the run quantification data
    if (any(is.element(colnames(rqall), "RUN"))) {
        rqall <- rqall[order(rqall$Protein, as.numeric(as.character(rqall$RUN))), ]
        rownames(rqall) <- NULL
    }
    
    #Mike: Below is for in-house verification occasionally
    #processedquant <- list(ProcessedData=work.NoImpute, RunlevelData=rqall, SummaryMethod=summaryMethod, ModelQC=rqmodelqc, PredictBySurvival=workpred, ImputedData=AbundanceAfterImpute)
    processedquant <- list(ProcessedData=work, 
                           RunlevelData=rqall, 
                           SummaryMethod=summaryMethod, 
                           ModelQC=rqmodelqc, 
                           PredictBySurvival=workpred)
    
    return(processedquant)
  
}

########################################################
# Manual function allowing foreach to return a list of multiple variables
resultsAsLists <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


########################################################
.runQuantification <- function(data, summaryMethod, 
                               equalFeatureVar, 
                               cutoffCensored, censoredInt, 
                               remove50missing, MBimpute, 
                               original_scale, logsum, 
                               featureSubset,
                               message.show) {
    
  ##Since the imputation has been done before feature selection, delete the columns of censoring indicator to avoid imputing the same intensity again   
  #if(featureSubset == "highQuality") {
  # data$cen <- NULL; data$pred <- NULL; data$INTENSITY <- 2^data$ABUNDANCE
  #} 
    
  ##If we want to impute again after the feature selection
  #if(featureSubset == "highQuality" & ImputeAgain==TRUE) {
  # data$ABUNDANCE <- data$ABUNDANCE.O  
  #}
   
    ## if there is 'remove' column, remove TRUE
    ## 2016. 08.29 should change column name for this remove variable. from feature selection??
    if( any(is.element(colnames(data), 'remove')) ) {
        data <- data[!data$remove, ]
    }

         
#   finalresult <- data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$RUN)),RUN=rep(c(levels(data$RUN)),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)

    # for saving predicting value for impute option
    predAbundance <- NULL
    
    ###################################
    ## method 1 : model based summarization
    if (summaryMethod == "linear"  & is.null(censoredInt)) {
        
        data <- data[!is.na(data$ABUNDANCE),]
        data$PROTEIN <- factor(data$PROTEIN)
        data$RUN <- factor(data$RUN)
    
        result <- NULL
        dataafterfit <- NULL
        
        for(i in 1: nlevels(data$PROTEIN)) {
        
            sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]

            sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
            sub$FEATURE <- factor(sub$FEATURE)  
            sub$RUN <- factor(sub$RUN)          
      
            temp <- data.frame(xtabs(~RUN, data=sub))

            sub.result <- data.frame(Protein=rep(unique(sub$PROTEIN),
                                         each=nlevels(sub$RUN)),
                                     RUN=rep(c(levels(sub$RUN)),1),
                                     LogIntensities=NA, 
                                     NumFeature=length(unique(sub$FEATURE)),
                                     NumPeaks=temp$Freq)
    
            singleFeature <- .checkSingleFeature(sub)
            singleSubject <- .checkSingleSubject(sub)
        
            ##### fit the model
            #if (message.show) {
                message(paste("Getting the summarization per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
            #}
            
            fit <- try(.fit.quantification.run(sub, singleFeature, singleSubject, equalFeatureVar), silent=TRUE)
             
            if (class(fit)=="try-error") {
                message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
                result <- rbind(result, sub.result)
          
                if (nrow(sub)!=0) {
                    sub$residuals <- NA
                    sub$fitted <- NA
                }
          
            } else {
          
                if (class(fit)=="lm") {
                    cf  <-  summary(fit)$coefficients
                }else{
                    cf  <-  fixef(fit)
                }
          
                # calculate sample quantification for all levels of sample
                a=1 
          
                for(j in 1:nlevels(sub$RUN)) {
              
                    contrast.matrix <- rep(0, nlevels(sub$RUN))
                    contrast.matrix[j] <- 1
            
                    contrast <- .make.contrast.run.quantification(fit,contrast.matrix,sub)
            
                    if (class(fit)=="lm") {
                        sub.result[a,3] <- .estimableFixedQuantification(cf,contrast)
                    } else {
                        sub.result[a,3] <- .estimableRandomQuantification(cf,contrast)
                    }
          
                    a=a+1
                }
          
                result <- rbind(result, sub.result)

                if (class(fit)=="lm") {  ### lm model
                    sub$residuals <- fit$residuals
                    sub$fitted <- fit$fitted.values
                } else {   ### lmer model
                    sub$residuals <- resid(fit)
                    sub$fitted <- fitted(fit)
                }
            
                dataafterfit <- rbind(dataafterfit,sub)

            }
        
        } ## end-loop for each protein  
    } ## for linear model summary
    
    ###################################
    ## Method 2 : Tukey Median Polish   
    if (summaryMethod == "TMP") {
        
        #data <- data[!is.na(data$ABUNDANCE),]
        data$PROTEIN <- factor(data$PROTEIN)
        data$RUN <- factor(data$RUN)
    
        result <- NULL
      
	pb <- progress_bar$new(format = "Progress [:bar] :percent eta :eta", total = nlevels(data$PROTEIN))
        
        for(i in 1: nlevels(data$PROTEIN)) {
            
            sub <- data[data$PROTEIN == levels(data$PROTEIN)[i], ]
                                  
            if (message.show) {
                message(paste("Getting the summarization by Tukey's median polish per subplot for protein ",
                              unique(sub$PROTEIN), "(", i," of ", length(unique(data$PROTEIN)), ")"))
            }
              
            sub$FEATURE <- factor(sub$FEATURE)  
              
            ### how to decide censored or not
            if ( MBimpute ) {
                ## 1. censored 
                if (censoredInt == "0") {
                      
                    sub[sub$censored == TRUE, 'ABUNDANCE'] <- 0
                    subtemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
                      
                }
                  
                ## 2. all censored missing
                if (censoredInt == "NA") {
                      
                    sub[sub$censored == TRUE, 'ABUNDANCE'] <- NA
                    subtemp <- sub[!is.na(sub$ABUNDANCE), ]
                      
                }
            } else {
                subtemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
            }
              
            ## if all measurements are NA,
            if ( nrow(sub) == (sum(is.na(sub$ABUNDANCE)) + sum(!is.na(sub$ABUNDANCE) & sub$ABUNDANCE == 0)) ) {
                message(paste("Can't summarize for ", unique(sub$PROTEIN), 
                              "(", i, " of ", length(unique(data$PROTEIN)),
                              ") because all measurements are NAs."))
                next()
            }
              
            countfeature <- xtabs(~FEATURE, subtemp)
            namefeature <- names(countfeature)[countfeature <= 1]
              
            if (length(namefeature) != 0) {
                sub <- sub[-which(sub$FEATURE %in% namefeature), ]
                  
                if (nrow(sub) == 0) {
                    message(paste("Can't summarize for ", unique(sub$PROTEIN), 
                            "(", i, " of ", length(unique(data$PROTEIN)), 
                            ") because features have one or fewer measurements across MS runs."))
                    next()
                      
                } else {
                    sub$FEATURE <- factor(sub$FEATURE)
                }
            }
             
            ## check one more time
            ## if all measurements are NA,
            if ( nrow(sub) == (sum(is.na(sub$ABUNDANCE)) + sum(!is.na(sub$ABUNDANCE) & sub$ABUNDANCE == 0)) ) {
                message(paste("After removing features which has only 1 measurement, Can't summarize for ",
                        unique(sub$PROTEIN), "(", i," of ", length(unique(data$PROTEIN)), 
                        ") because all measurements are NAs."))
                next()
            }
              
            ## remove run which has no measurement at all 
            ## remove features which are completely NAs
            if ( MBimpute ) {
                ## 1. censored 
                if (censoredInt == "0") {
                    subtemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
                }
                  
                ## 2. all censored missing
                if (censoredInt == "NA") {
                    subtemp <- sub[!is.na(sub$ABUNDANCE), ]
                }  
            } else {
                subtemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
            }
              
            count <- aggregate(ABUNDANCE ~ RUN, data=subtemp, length)
            norun <- setdiff(unique(data$RUN), count$RUN)
              
            if (length(norun) != 0 & length(intersect(norun, as.character(unique(sub$RUN))))) { 
                # removed NA rows already, if there is no overlapped run, error
                sub <- sub[-which(sub$RUN %in% norun), ]
                sub$RUN <- factor(sub$RUN)
            }
              
            if (remove50missing) {
                # count # feature per run
                if (!is.null(censoredInt)) {
                    if (censoredInt == "NA") {
                        subtemp <- sub[!is.na(sub$INTENSITY), ]
                    }
                      
                    if (censoredInt == "0") {
                        subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
                    }
                }
                  
                numFea <- xtabs(~RUN, subtemp) ## RUN or run.label?
                numFea <- numFea/length(unique(subtemp$FEATURE))
                numFea <- numFea <= 0.5
                removerunid <- names(numFea)[numFea]
                  
                ## if all measurements are NA,
                if (length(removerunid)==length(numFea)) {
                    message(paste("Can't summarize for ",unique(sub$PROTEIN), 
                            "(", i, " of ", length(unique(data$PROTEIN)),
                            ") because all runs have more than 50% NAs and are removed with the option, remove50missing=TRUE."))
                    next()
                }
                  
            }
              
            ## check whether we need to impute or not.
            if (sum(sub$censored) > 0) {
                  
                ## 2. put minimum in feature level to NA
                if (cutoffCensored == "minFeature") {
                    if (censoredInt == "NA") {
                        cut <- aggregate(ABUNDANCE ~ FEATURE, data=sub, function(x) min(x, na.rm=TRUE))
                        ## cutoff for each feature is little less than minimum abundance in a run.
                        cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
                          
                        ## remove runs which has more than 50% missing values
                        if (remove50missing) {
                            if (length(removerunid) != 0) {
                                sub <- sub[-which(sub$RUN %in% removerunid), ]
                                sub$RUN <- factor(sub$RUN)
                            }
                        }
                          
                        for(j in 1:length(unique(cut$FEATURE))) {
                            sub[is.na(sub$ABUNDANCE) & sub$FEATURE == cut$FEATURE[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                        }
                    }
                      
                    if (censoredInt == "0") {
                        subtemptemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
                        cut <- aggregate(ABUNDANCE ~ FEATURE, data=subtemptemp, FUN=min)
                        ## cutoff for each feature is little less than minimum abundance in a run.
                        cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
                          
                        ## remove runs which has more than 50% missing values
                        if (remove50missing) {
                            if (length(removerunid) != 0) {
                                sub <- sub[-which(sub$RUN %in% removerunid), ]
                                sub$RUN <- factor(sub$RUN)
                            }
                        }
                        
                        for(j in 1:length(unique(cut$FEATURE))) {
                            sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE == 0  & 
                                    sub$FEATURE == cut$FEATURE[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                        }
                    }
                }
                  
                ## 3. put minimum in RUN to NA
                if (cutoffCensored == "minRun") {
                      
                    ## remove runs which has more than 50% missing values
                    if (remove50missing) {
                        if (length(removerunid) != 0) {
                            sub <- sub[-which(sub$RUN %in% removerunid), ]
                            sub$RUN <- factor(sub$RUN)
                        }
                    }
                      
                    if (censoredInt == "NA") {
                        cut <- aggregate(ABUNDANCE ~ RUN, data=sub, function(x) min(x, na.rm=TRUE))
                        ## cutoff for each Run is little less than minimum abundance in a run.
                        cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
                          
                        for(j in 1:length(unique(cut$RUN))) {
                            sub[is.na(sub$ABUNDANCE) & 
                                    sub$RUN == cut$RUN[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                        }
                    }
                      
                    if (censoredInt == "0") {
                        subtemptemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
                        cut <- aggregate(ABUNDANCE ~ RUN, data=subtemptemp, FUN=min)
                        cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
                          
                        for(j in 1:length(unique(cut$RUN))) {
                            sub[!is.na(sub$ABUNDANCE) & 
                                    sub$ABUNDANCE == 0 & 
                                    sub$RUN == cut$RUN[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                        }
                    }
                }
                  
                ## 20150829 : 4. put minimum RUN and FEATURE
                if (cutoffCensored == "minFeatureNRun") {
                    if (censoredInt == "NA") {
                          
                        ## cutoff for each feature is little less than minimum abundance in a run.
                          
                        cut.fea <- aggregate(ABUNDANCE ~ FEATURE, data=sub, function(x) min(x, na.rm=TRUE))
                        cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
                          
                        ## remove runs which has more than 50% missing values
                        ## before removing, need to contribute min feature calculation
                        if (remove50missing) {
                            if (length(removerunid) != 0) {
                                sub <- sub[-which(sub$RUN %in% removerunid), ]
                                sub$RUN <- factor(sub$RUN)
                            }
                        }
                          
                        ## cutoff for each Run is little less than minimum abundance in a run.
                          
                        cut.run <- aggregate(ABUNDANCE ~ RUN, data=sub, function(x) min(x, na.rm=TRUE))
                        cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
                          
                        if (length(unique(cut.fea$FEATURE)) > 1) {
                            for(j in 1:length(unique(cut.fea$FEATURE))) {
                                for(k in 1:length(unique(cut.run$RUN))) {
                                    # get smaller value for min Run and min Feature
                                    finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
                                      
                                    sub[is.na(sub$ABUNDANCE) & 
                                            sub$FEATURE == cut.fea$FEATURE[j] & 
                                            sub$RUN == cut.run$RUN[k], "ABUNDANCE"] <- finalcut
                                }
                            }
                        }
                        # if single feature, not impute
                    }
                      
                    if (censoredInt == "0") {
                        subtemptemp <- sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE != 0, ]
                          
                        cut.fea <- aggregate(ABUNDANCE ~ FEATURE, data=subtemptemp, FUN=min)
                        cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
                          
                        ## remove runs which has more than 50% missing values
                        ## before removing, need to contribute min feature calculation
                        if (remove50missing) {
                            if (length(removerunid) != 0) {
                                sub <- sub[-which(sub$RUN %in% removerunid), ]
                                sub$RUN <- factor(sub$RUN)
                            }
                        }
                          
                        cut.run <- aggregate(ABUNDANCE~RUN, data=subtemptemp, FUN=min)
                        cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
                          
                        if (length(unique(cut.fea$FEATURE)) > 1) {
                            for(j in 1:length(unique(cut.fea$FEATURE))) {
                                for(k in 1:length(unique(cut.run$RUN))) {
                                    # get smaller value for min Run and min Feature
                                    finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                                      
                                    sub[!is.na(sub$ABUNDANCE) & 
                                            sub$ABUNDANCE == 0 & 
                                            sub$FEATURE == cut.fea$FEATURE[j] & 
                                            sub$RUN == cut.run$RUN[k], "ABUNDANCE"] <- finalcut
                                }
                            }
                        } else { # single feature
                              
                            sub[!is.na(sub$ABUNDANCE) & sub$ABUNDANCE == 0, "ABUNDANCE"] <- cut.fea$ABUNDANCE
                              
                        }
                    }
                }
                  
                if (MBimpute) {
                      
                    if (nrow(sub[sub$censored, ]) > 0) {
                        ## impute by survival model
                        subtemp <- sub[!is.na(sub$ABUNDANCE),]
                        countdf <- nrow(subtemp) < (length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
                          
                        set.seed(100)
                        ### fit the model
                        if (length(unique(sub$FEATURE)) == 1) {
                            fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN, 
                                                           data=sub, dist='gaussian')
                        } else {
                            if (countdf) {
                                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN, 
                                                             data=sub, dist='gaussian')
                            } else {
                                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN, 
                                                               data=sub, dist='gaussian')
                            }
                        }
                          
                        # get predicted value from survival
                        sub <- data.frame(sub, pred=predict(fittest, newdata=sub, type="response"))
                          
                        # the replace censored value with predicted value
                        sub[sub$censored, "ABUNDANCE"] <- sub[sub$censored, "pred"] 
                          
                        # save predicted value
                          # predAbundance <- c(predAbundance,predict(fittest, newdata=sub, type="response"))
                          #predAbundance <- c(predict(fittest, newdata=sub, type="response"))
                    }
                }   
            }
              
            ## then, finally remove NA in abundance
            sub <- sub[!is.na(sub$ABUNDANCE), ]
              
            if (nlevels(sub$FEATURE) > 1) { ## for more than 1 features
                  
                data_w <- reshape2::dcast(RUN ~ FEATURE, data=sub, value.var='ABUNDANCE', keep=TRUE)
                rownames(data_w) <- data_w$RUN
                data_w <- data_w[, -1]
                data_w[data_w == 1] <- NA
                  
                if (!original_scale) {
                      
                    meddata  <-  medpolish(data_w,na.rm=TRUE, trace.iter = FALSE)
                    tmpresult <- meddata$overall + meddata$row
                      
                      ## if fractionated sample, need to get per sample run
                      ## ?? if there are technical replicates, how to match sample and MS run for different fractionation??
                      
                      #if( length(unique(sub$METHOD)) > 1 ) {
                      #  runinfo <- unique(sub[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN", "METHOD")])
                      #  runinfo$uniquesub <- paste(runinfo$GROUP_ORIGINAL, runinfo$SUBJECT_ORIGINAL, sep="_")
                      #}
                      
                } else { # original_scale
                    data_w <- 2^(data_w)
                     
                    meddata  <-  medpolish(data_w,na.rm=TRUE, trace.iter = FALSE)
                    tmpresult <- meddata$overall + meddata$row
                      
                    tmpresult <- log2(tmpresult)
                }
                  
                # count # feature per run
                if (!is.null(censoredInt)) {
                    if (censoredInt == "NA") {
                        subtemp <- sub[!is.na(sub$INTENSITY), ]
                        subtempimpute <- sub[is.na(sub$INTENSITY), ]
                        subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE), ]
                    }
                      
                    if (censoredInt == "0") {
                        subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
                        subtempimpute <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0, ]
                        subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE) & subtempimpute$ABUNDANCE != 0, ]
                    }
                      
                    subtemp$RUN <- factor(subtemp$RUN, levels = rownames(data_w))
                    numFea <- xtabs(~RUN, subtemp)
                    numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
                    numFeaTF <- numFeaPercentage >= 0.5
                      
                    subtempimpute$RUN <- factor(subtempimpute$RUN, levels = rownames(data_w))
                    numimpute <- xtabs(~RUN, subtempimpute)
                      
                    sub.result <- data.frame(Protein = unique(sub$PROTEIN), 
                                             LogIntensities = tmpresult, 
                                             RUN = names(tmpresult), 
                                             NumMeasuredFeature = as.vector(numFea), 
                                             MissingPercentage = as.vector(numFeaPercentage), 
                                             more50missing = numFeaTF, 
                                             NumImputedFeature = as.vector(numimpute))
                      
                } else {
                    subtemp <- sub[!is.na(sub$INTENSITY), ]
                      
                    subtemp$RUN <- factor(subtemp$RUN, levels =rownames(data_w))
                    numFea <- xtabs(~RUN, subtemp)
                    numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
                    numFeaTF <- numFeaPercentage >= 0.5
                      
                    sub.result <- data.frame(Protein=unique(sub$PROTEIN), 
                                             LogIntensities=tmpresult, 
                                             RUN=names(tmpresult), 
                                             NumMeasuredFeature = as.vector(numFea), 
                                             MissingPercentage=as.vector(numFeaPercentage), 
                                             more50missing=numFeaTF)
                      
                }
                  
                result <- rbind(result, sub.result)
            } else { ## single feature
                  
                ## single feature, use original values
                  
                subtemp <- sub[!is.na(sub$ABUNDANCE),]
                  
                if (!is.null(censoredInt)) {
                    if (censoredInt == "NA") {
                        subtempcount <- sub[!is.na(sub$INTENSITY), ]
                        subtempimpute <- sub[is.na(sub$INTENSITY), ]
                        subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE), ]
                    }
                      
                    if (censoredInt == "0") {
                        subtempcount <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
                        subtempimpute <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0, ]
                        subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE) & subtempimpute$ABUNDANCE != 0, ]
                    }
                      
                    numFea <- xtabs(~RUN, subtempcount)
                    numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
                    numFeaTF <- numFeaPercentage >= 0.5
                      
                    numimpute <- xtabs(~RUN, subtempimpute)
                      
                    sub.result <- data.frame(Protein=subtemp$PROTEIN,
                                             LogIntensities=subtemp$ABUNDANCE, 
                                             RUN=subtemp$RUN, 
                                             NumMeasuredFeature = as.vector(numFea), 
                                             MissingPercentage=as.vector(numFeaPercentage), 
                                             more50missing=numFeaTF, 
                                             NumImputedFeature = as.vector(numimpute))
                      
                } else {
                    subtempcount <- subtemp
                      
                    numFea <- xtabs(~RUN, subtempcount)
                    numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
                    numFeaTF <- numFeaPercentage >= 0.5
                      
                    sub.result <- data.frame(Protein=subtemp$PROTEIN,
                                             LogIntensities=subtemp$ABUNDANCE, 
                                             RUN=subtemp$RUN, 
                                             NumMeasuredFeature = as.vector(numFea), 
                                             MissingPercentage=as.vector(numFeaPercentage), 
                                             more50missing=numFeaTF)
                      
                }
                  
                result <- rbind(result, sub.result)
            }
            
            ## progress
            pb$tick()
            
        }  ## loop for proteins
        pb$terminate()
        
        dataafterfit <- NULL
    }
        
    ###################################
    ## Method 3 : log sum   
    ## retired on Aug 2 2016
    
    ###################################
    ## method 4 : survival model for censored missing values
    if (summaryMethod == "linear" & !is.null(censoredInt)) {
            
        #data <- data[!is.na(data$ABUNDANCE),]
        data$PROTEIN <- factor(data$PROTEIN)
        data$RUN <- factor(data$RUN)
        
        result <- NULL

        for(i in 1:length(unique(data$PROTEIN))) {

            sub <- data[data$PROTEIN == unique(data$PROTEIN)[i], ]
            
            if (message.show) {
                message(paste("Getting the summarization for censored missing values per subplot for protein ",
                              unique(sub$PROTEIN), "(", i, " of ", length(unique(data$PROTEIN)), ")"))
            }
            
            sub$FEATURE <- factor(sub$FEATURE)

            ## if all measurements are NA,
            if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
                message(paste("Can't summarize for ", unique(sub$PROTEIN), 
                              "(", i, " of ", length(unique(data$PROTEIN)),
                              ") because all measurements are NAs."))
                next()
            }
            
            ## remove run which has no measurement at all
            subtemp <- sub[!is.na(sub$INTENSITY), ]
            count <- aggregate(ABUNDANCE~RUN, data=subtemp, length)
            norun <- setdiff(unique(data$RUN), ount$RUN)
            
            if (length(norun) != 0 & length(intersect(norun, as.character(unique(sub$RUN)))) != 0) { 
                # removed NA rows already, if there is no overlapped run, error
                sub <- sub[-which(sub$RUN %in% norun), ]
                sub$RUN <- factor(sub$RUN)
            }
            
            if (length(unique(sub$RUN)) == 1) {
            
                message(paste("* Only 1 MS run in ", levels(data$PROTEIN)[i], 
                              " has measurement. Can't summarize with censored intensities."))
            
                next()
            }   
                    
            
            ## remove features which are (completely NAs or zero) 
            subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
            countfeature <- xtabs(~FEATURE, subtemp)
            namefeature <- names(countfeature)[countfeature == 0]
            
            if (length(namefeature) != 0) {
                sub <- sub[-which(sub$FEATURE %in% namefeature), ]
                sub$FEATURE <- factor(sub$FEATURE)
            }
            
            if (nrow(sub) == 0) {
            
                message(paste("* All measurements are NAs or only one measurement per feature in ",
                              levels(data$PROTEIN)[i], ". Can't summarize with censored intensities."))
            
                next()
            }   

            ##### how to decide censored or not
            ## 1. censored 
            if (censoredInt == "0") {
                sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY == 0, 0, 1)
            }
            
            ### 2. all censored missing
            if (censoredInt == "NA") {
                sub$cen <- ifelse(is.na(sub$INTENSITY), 0, 1)
            }

            ## 2. put minimum in feature level to NA
            if (cutoffCensored == "minFeature") {
                if (censoredInt == "NA") {
                    cut <- aggregate(ABUNDANCE ~ FEATURE, data=sub, function(x) min(x, na.rm=TRUE))
                    ## cutoff for each Run is little less than minimum abundance in a run.
                    cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

                    for(j in 1:length(unique(cut$FEATURE))) {
                        sub[is.na(sub$INTENSITY) & sub$FEATURE == cut$FEATURE[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                    }
                }
                
                if (censoredInt == "0") {
                    subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
                    cut <- aggregate(ABUNDANCE ~ FEATURE, data=subtemptemp, FUN=min)
                    ## cutoff for each Run is little less than minimum abundance in a run.
                    cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

                    for(j in 1:length(unique(cut$FEATURE))) {
                        sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0 & 
                                sub$FEATURE == cut$FEATURE[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                    }
                }
            }
            
            ## 3. put minimum in RUN to NA
            if (cutoffCensored == "minRun") {
                if (censoredInt == "NA") {
                    cut <- aggregate(ABUNDANCE~RUN, data=sub, function(x) min(x, na.rm=TRUE))
                    ## cutoff for each Run is little less than minimum abundance in a run.
                    cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

                    for(j in 1:length(unique(cut$RUN))) {
                        sub[is.na(sub$INTENSITY) & sub$RUN == cut$RUN[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                    }
                }
                
                if (censoredInt == "0") {
                    subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
                    cut <- aggregate(ABUNDANCE~RUN, data=subtemptemp, FUN=min)
                    
                    ## cutoff for each Run is little less than minimum abundance in a run.
                    cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

                    for(j in 1:length(unique(cut$RUN))) {
                        sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0 & sub$RUN==cut$RUN[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
                    }
                }
            }   
            
            ## 20150829 : 4. put minimum RUN and FEATURE
            if (cutoffCensored == "minFeatureNRun") {
                if (censoredInt == "NA") {
                    
                    ## cutoff for each feature is little less than minimum abundance in a run.
                    cut.fea <- aggregate(ABUNDANCE ~ FEATURE, data=sub, function(x) min(x, na.rm=TRUE))
                    cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
                                            
                    ## cutoff for each Run is little less than minimum abundance in a run.

                    cut.run <- aggregate(ABUNDANCE ~ RUN, data=sub, function(x) min(x, na.rm=TRUE))
                    cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
                    
                    if (length(unique(sub$FEATURE)) > 1) {
                        for(j in 1:length(unique(sub$FEATURE))) {
                            for(k in 1:length(unique(sub$RUN))) {
                                # get smaller value for min Run and min Feature
                                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                            
                                sub[is.na(sub$INTENSITY) & sub$FEATURE == cut.fea$FEATURE[j] & 
                                        sub$RUN == cut.run$RUN[k], "ABUNDANCE"] <- finalcut
                            }
                        }
                    }
                        # if single feature, not impute
                }
                
                if (censoredInt == "0") {
                    subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]

                    cut.fea <- aggregate(ABUNDANCE ~ FEATURE, data=subtemptemp, FUN=min)
                    cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
                    
                    cut.run <- aggregate(ABUNDANCE ~ RUN, data=subtemptemp, FUN=min)
                    cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE

                    if (length(unique(sub$FEATURE)) > 1) {
                        for(j in 1:length(unique(sub$FEATURE))) {
                            for(k in 1:length(unique(sub$RUN))) {
                                # get smaller value for min Run and min Feature
                                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                            
                                sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0 & 
                                        sub$FEATURE == cut.fea$FEATURE[j] & sub$RUN == cut.run$RUN[k], "ABUNDANCE"] <- finalcut
                            }
                        }
                    } else { # single feature
                
                        sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0, "ABUNDANCE"] <- cut.fea$ABUNDANCE
                        
                    }
                }
            }
                
            
            ## when number of measurement is less than df, error for fitting
            subtemp <- sub[!is.na(sub$ABUNDANCE), ]
            countdf <- nrow(subtemp) < (length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
            
            ### fit the model
            if (length(unique(sub$FEATURE)) == 1) {
                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN,
                                             data=sub, dist='gaussian')
            } else {
                if (countdf) {
                    fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN,
                                                 data=sub, dist='gaussian')
                } else {
                    fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN,
                                                 data=sub, dist='gaussian')
                }
            }
            
            sub.result <- data.frame(Protein=unique(sub$PROTEIN),
                                     RUN=rep(c(levels(sub$RUN)), 1),
                                     LogIntensities=NA)

            # get the parameters
            cf <- summary(fittest)$coefficients

            # calculate sample quantification for all levels of sample
            a <- 1  
      
            for(j in 1:nlevels(sub$RUN)) {
                contrast.matrix <- rep(0, nlevels(sub$RUN))
                contrast.matrix[j] <- 1
                        contrast <- .make.contrast.run.quantification.Survival(fittest, contrast.matrix,sub)

                 sub.result[a, 3] <- .estimableFixedQuantificationSurvival(cf, contrast)
                 a <- a+1
            }

            result <- rbind(result, sub.result)
        }

        datamat <- reshape2::dcast( Protein ~ RUN, data=result, value.var='LogIntensities', keep=TRUE)
        datamat <- melt(datamat, id.vars=c('Protein'))
        colnames(datamat) <- c('Protein','RUN','LogIntensities')
        result <- datamat
    }
    dataafterfit <- NULL    
    
    ###################################
    ## final result
    finalout <- list(rqdata=result, ModelQC=dataafterfit, PredictedBySurvival=predAbundance)
    return(finalout)
}



##########################################################################################
## updated v3
.fit.quantification.run <- function(sub, singleFeature, singleSubject, equalFeatureVar) {
    
    ## for single Feature, original value is the run quantification
    if (singleFeature) {
        fit.full <- lm(ABUNDANCE ~ RUN , data = sub)
    }else{
        fit.full <- lm(ABUNDANCE ~ FEATURE + RUN , data = sub)
    }
    
    
    ## make equal variance for feature : need to be updated
    if (!equalFeatureVar) {
       fit.full <- .iter.wls.fit.model(data=sub, fit=fit.full, nrepeats=1)
    }
    
    return(fit.full)
}
