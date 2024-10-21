### Taxonomy
allSubNames <- c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus",
                 "Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
numberLetters <- length(sort(unique(unlist(strsplit(allSubNames,""," ")))))

path_structure_dino <- list(
  Pisanosaurus = list(p1 = c(0,0,0)),
  Rinchenia = list(p1=c(0,0,1)),
  Edmontosaurus = list(p1=c(0,1,0)),
  Parasaurolophus = list(p1=c(0,1,1)),
  Rhamphorhynchinae = list(p1=c(1,0,0),p2=c(1,0,1)),
  Pteranodontia = list(p1=c(1,1,0),p2=c(1,1,1))
)
path_structure_ptero <- list(
  Ornithoscelida = list(p1=c(0,0,0),p2=c(0,0,1)),
  Saurolophidae = list(p1=c(0,1,0),p2=c(0,1,1)),
  Dorygnathus = list(p1=c(1,0,0)),
  Scaphognathus = list(p1=c(1,0,1)),
  Istiodactylus = list(p1=c(1,1,0)),
  Pteranodon = list(p1=c(1,1,1))
)

cat_neighbors <- list(
  Pisanosaurus = list(
    correct_on_L3 = c("Pisanosaurus"),
    correct_on_L2 = c("Pisanosaurus","Rinchenia"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Rinchenia = list(
    correct_on_L3 = c("Rinchenia"),
    correct_on_L2 = c("Pisanosaurus","Rinchenia"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Edmontosaurus = list(
    correct_on_L3 = c("Edmontosaurus"),
    correct_on_L2 = c("Edmontosaurus","Parasaurolophus"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Parasaurolophus = list(
    correct_on_L3 = c("Parasaurolophus"),
    correct_on_L2 = c("Edmontosaurus","Parasaurolophus"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Dorygnathus = list(
    correct_on_L3 = c("Dorygnathus"),
    correct_on_L2 = c("Dorygnathus","Scaphognathus"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Scaphognathus = list(
    correct_on_L3 = c("Scaphognathus"),
    correct_on_L2 = c("Dorygnathus","Scaphognathus"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Istiodactylus = list(
    correct_on_L3 = c("Istiodactylus"),
    correct_on_L2 = c("Istiodactylus","Pteranodon"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Pteranodon = list(
    correct_on_L3 = c("Pteranodon"),
    correct_on_L2 = c("Istiodactylus","Pteranodon"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  )
)

correctLabels_perLevel_dino <- list(
  Pisanosaurus = list(
    correct_on_L3 = c("Pisanosaurus"),
    correct_on_L2 = c("Pisanosaurus","Rinchenia"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Rinchenia = list(
    correct_on_L3 = c("Rinchenia"),
    correct_on_L2 = c("Pisanosaurus","Rinchenia"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Edmontosaurus = list(
    correct_on_L3 = c("Edmontosaurus"),
    correct_on_L2 = c("Edmontosaurus","Parasaurolophus"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Parasaurolophus = list(
    correct_on_L3 = c("Parasaurolophus"),
    correct_on_L2 = c("Edmontosaurus","Parasaurolophus"),
    correct_on_L1 = c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus")
  ),
  Rhamphorhynchinae = list(
    correct_on_L3 = c("Rhamphorhynchinae"),
    correct_on_L2 = c("Rhamphorhynchinae"),
    correct_on_L1 = c("Rhamphorhynchinae","Pteranodontia")
  ),
  Pteranodontia = list(
    correct_on_L3 = c("Pteranodontia"),
    correct_on_L2 = c("Pteranodontia"),
    correct_on_L1 = c("Rhamphorhynchinae","Pteranodontia")
  )
)

correctLabels_perLevel_ptero <- list(
  Ornithoscelida = list(
    correct_on_L3 = c("Ornithoscelida"),
    correct_on_L2 = c("Ornithoscelida"),
    correct_on_L1 = c("Ornithoscelida","Saurolophidae")
  ),
  Saurolophidae = list(
    correct_on_L3 = c("Saurolophidae"),
    correct_on_L2 = c("Saurolophidae"),
    correct_on_L1 = c("Ornithoscelida","Saurolophidae")
  ),
  Dorygnathus = list(
    correct_on_L3 = c("Dorygnathus"),
    correct_on_L2 = c("Dorygnathus","Scaphognathus"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Scaphognathus = list(
    correct_on_L3 = c("Scaphognathus"),
    correct_on_L2 = c("Dorygnathus","Scaphognathus"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Istiodactylus = list(
    correct_on_L3 = c("Istiodactylus"),
    correct_on_L2 = c("Istiodactylus","Pteranodon"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  ),
  Pteranodon = list(
    correct_on_L3 = c("Pteranodon"),
    correct_on_L2 = c("Istiodactylus","Pteranodon"),
    correct_on_L1 = c("Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  )
)

code_fx <- function(file){
	data <- read.table(file,sep=",",fill=T,head = T)
	code <- c(as.matrix(data[3,"responses"]))
	code <- unlist(strsplit(code,":"))[2]
	code <- toupper(gsub("\\}|\\}|\"","",code))
	return(code)
}

fxOspan <- function(file){
	
	data <- read.table(file,sep=",",fill=T,head = T)
	
	code <- c(as.matrix(data[3,"responses"]))
	code <- unlist(strsplit(code,":"))[2]
	code <- toupper(gsub("\\}|\\}|\"","",code))
	
	summaryData <- data[which(data[,"summary"]=="summary_realTask"),]
	
	solve_accuracy <- mean(summaryData$solve_accuracy)
	
	nCorrect <- sum(summaryData[,"nLettersCorrect"])
	nPresented <- sum(summaryData[,"nLettersPresented"])
	propCorrect <- nCorrect/nPresented
	
	toBeReturned <- data.frame(code, solve_accuracy, nCorrect, nPresented, propCorrect)
	
	return (toBeReturned)
	
}

taxoFx <- function(file){
  data <- read.table(file,sep=",",fill=T,head = T)
  code <- c(as.matrix(data[3,"responses"]))
  code <- unlist(strsplit(code,":"))[2]
  code <- toupper(gsub("\\}|\\}|\"","",code))
  data_judgement <- data[which(data$summary_judgment =="summary_judgment"),]
  presentedLabels <- as.matrix(data_judgement$presentedNamesLabels)[1]
  presentedLabels <- unlist(strsplit(presentedLabels,","))
  if(is.element("Edmontosaurus", presentedLabels)){taxo.type="taxo.dino"}else{taxo.type="taxo.ptero"}
  return(data.frame(file,code,taxo.type))
}

"fxDinoNimi <- function(file){
	
	#data <- read.table(listedFilesDinoNimi[76],sep=',',fill=T,head = T)
	data <- read.table(file,sep=',',fill=T,head = T)
	
	code <- c(as.matrix(data[3,'responses']))
	code <- unlist(strsplit(code,':'))[2]
	code <- toupper(gsub('\\}|\\}|\'','',code))
	
	data_judgement <- data[which(data[,'summary_judgment'] =='summary_judgment'),]
	
	mean.rt <- mean(data_judgement$rt[1])
	
	judgments_per_phase <- strsplit(as.matrix(data_judgement[,'judgments']),',')
	judgments_per_phase <- lapply(judgments_per_phase,function(i){as.numeric(i)+1})
	judgmentProps_per_phase <- lapply(judgments_per_phase,function(i){
		i_freqs <- prop.table(table(i))
		i_freqs_vec <- rep(0,5)
		i_freqs_vec[as.numeric(names(i_freqs))] <- i_freqs
		return(i_freqs_vec)
		}
		
		)
	judgmentProps_per_phase <- do.call(rbind, judgmentProps_per_phase)
	nRounds <- nrow(judgmentProps_per_phase)
	rownames(judgmentProps_per_phase) <- paste('Judgment-Round',c(1: nRounds),sep='')
	colnames(judgmentProps_per_phase) <- c('Sure yes','Prob. yes','Dunno','Prob. no','Sure no')
	
	data_categorization <- data[which(data[,'summary_categorization'] =='summary_categorization'),]
	RightWrong_per_phase <- table(data_categorization[,'genPhase'],data_categorization[,'correctCategorization'])
	if(length(which(data_categorization[,'correctCategorization']==0))==0){
		RightWrong_per_phase <- cbind(rep(0,nrow(RightWrong_per_phase)),RightWrong_per_phase)
		colnames(RightWrong_per_phase) <- c(0,1)
	}
	nSelected_per_Phase <- rowSums(RightWrong_per_phase)
	RW_per_phase_props <- prop.table(RightWrong_per_phase,1)
	
	consecStimuliSelected <- gsub('\\d','',data_categorization[,'stimulus'])
	for(i in 1:length(consecStimuliSelected)){
		if(i==1){nSelected <- NULL}
		nUnique_i <- length(unique(consecStimuliSelected[1:i]))
		nSelected <- c(nSelected, nUnique_i)
	}
	data_categorization <- cbind(data_categorization, nSelected)
	nPhases_select <- max(data_categorization[,'genPhase'])
	nNewSelected_per_phase <- sapply(c(1: nPhases_select),function(i){max(data_categorization[data_categorization[,'genPhase']==i,'nSelected'])})
	
	nUniqueSelected_per_phase <- sapply(c(1: nPhases_select), data_categorization= data_categorization,function(i, data_categorization){
		iData <- data_categorization <- data_categorization[data_categorization[,'genPhase']==i,]
		iData_stimuli <- gsub('\\d','', iData[,'stimulus'])
		return(length(unique(iData_stimuli)))
		})
	
	
	RW_per_phase_props <- cbind(nSelected_per_Phase, nUniqueSelected_per_phase ,nNewSelected_per_phase, RightWrong_per_phase ,RW_per_phase_props)
	colnames(RW_per_phase_props) <- c('nExp','nUnique','Cumul','nWrong','nRight','pWrong','pRight')
	nRounds <- nrow(RW_per_phase_props)
	rownames(RW_per_phase_props) <- paste('Selec-R',c(1: nRounds),sep='')
	
	toBeReturned <- list(judgmentProps_per_phase, RW_per_phase_props, mean.rt)
	names(toBeReturned) <- paste(code,'-',c('JudgmentBehavior','SelectionBehavior','Mean.RT'),sep='')
	
	return (toBeReturned)
	
}"

{fx_dinoNimi_judgment <- function(file, nRatingSteps, oSpan_Scores_df,
                                 correctLabels_perLevel_dino, correctLabels_perLevel_ptero){
  
  data <- read.table(file,sep=",",fill=T,head = T)
  
  code <- c(as.matrix(data[3,"responses"]))
  code <- unlist(strsplit(code,":"))[2]
  code <- toupper(gsub("\\}|\\}|\"","",code))
  
  oSpan.score <- oSpan_Scores_df$nCorrect[oSpan_Scores_df$code==code]
  WMC.group <- oSpan_Scores_df$WMC.group[oSpan_Scores_df$code==code]
  
  data_judgement <- data[which(data$summary_judgment =="summary_judgment"),]
  
  # Determining which exemplars were learned on category vs. item level
  presentedLabels <- as.matrix(data_judgement$presentedNamesLabels)
  presentedLabels <- strsplit(presentedLabels,",")
  if(is.element("Edmontosaurus", unlist(presentedLabels))){taxo.type="taxo.dino"}else{taxo.type="taxo.ptero"}
  # taxo-specific labels  
  if(taxo.type=="taxo.dino"){
    #correctLabels_per_level[[h]][["correct_on_L2"]]
    correctLabels_per_level <- correctLabels_perLevel_dino
  }else{
    correctLabels_per_level <- correctLabels_perLevel_ptero
  }
  
  
  timeStamps <- data_judgement[,"time_elapsed"]
  matchesNomatches <- as.matrix(data_judgement[,"matchNomatch"])
  matchesNomatches <- t(apply(matchesNomatches,1,function(i){as.numeric(unlist(strsplit(i,",")))}))
  providedResponses <- as.matrix(data_judgement[,"judgments"])
  providedResponses <- t(apply(providedResponses,1,function(i){as.numeric(unlist(strsplit(i,",")))}))
  
  presentedStimuli <- strsplit(as.matrix(data_judgement$presentedPicturesLabels),",")
  
  nJudgePhases <- nrow(providedResponses)
  nQuestions <- ncol(providedResponses)
  
  # for loop to recode responses: whether the responses of e.g. "sure yes" or "sure no" are correct depends on
  # whether the presented picture and label have matched. A correct response provided with high confidence (e.g.
  # "sure yes" when picture and label have matched, is recoded as 5, a correct response with low confidence
  # as 4, and so on)
  
  # js: 0=sure yes, 4 = sure no
  
  for(t in 1: nJudgePhases){  # each row corresponds to one judgment phase
    if(t==1){
      providedResponses_recoded <- list(L1=list(),L2=list(),L3=list())
      sdt.mat <- matrix(0,nrow=nJudgePhases,ncol=5)
      colnames(sdt.mat) <- c("Hit","FA","CR","Miss","Dunno")
      rownames(sdt.mat) <- paste("Phase_",c(1:nJudgePhases),sep="")
      sdt.mat <- list(L1=sdt.mat,L2=sdt.mat,L3=sdt.mat)
    }
    t_matchesNomatches <- matchesNomatches[t,]
    t_providedResponses <- providedResponses[t,]
    t_presentedStimuli <- presentedStimuli[[t]]
    t_labels <- presentedLabels[[t]]
    #correctLabels_per_level=correctLabels_per_level;w=t_matchesNomatches; x=t_labels;y=t_presentedStimuli;z=t_providedResponses
    recodeFx <- function(correctLabels_per_level, w, x, y, z){
      for(i in 1: nQuestions){
        if(i==1){
          providedResponses_recoded_t <- list(L1=list(),L2=list(),L3=list())
          sdt.mat <- matrix(0,nrow=nJudgePhases,ncol=5)
          colnames(sdt.mat) <- c("Hit","FA","CR","Miss","Dunno")
          rownames(sdt.mat) <- paste("Phase_",c(1:nJudgePhases),sep="")
          sdt.mat_t <- list(L1=sdt.mat,L2=sdt.mat,L3=sdt.mat)
        }
        iLabel <- x[i]
        iPicture <- y[i]
        iResponse <- z[i]
        iMatch <- w[i]
        # Whether response is correct depends on the level of observation. However, if response is correct on 
        # Level i, then it is also correct on Levels >i
        choice_correct <- FALSE
        for(L_x in 1:3){
          if(choice_correct==TRUE){
            iResponse_recoded <- iResponse_recoded
            sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1
          }else{
            allowedLabelChoices <- correctLabels_per_level[[iPicture]][[L_x]]
            if(is.element(iLabel,allowedLabelChoices)){
              L_x.Match <- 1
            }else{L_x.Match <- 0}
            
            if(L_x.Match==1){
              recodeVec <- nRatingSteps:1
            }else{recodeVec <- 1:nRatingSteps}
            
            iResponse_recoded <- recodeVec[iResponse+1]
            if(iResponse_recoded==3){
              iResponse_sdt <- "Dunno"
              sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1
            }else{
              if(iResponse_recoded<=2){
                if(L_x.Match==0){
                  iResponse_sdt <- "FA"
                  sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1
                }else{
                  iResponse_sdt <- "Miss"
                  sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1}
              }
              if(iResponse_recoded>=4){
                choice_correct <- TRUE
                if(L_x.Match==0){
                  iResponse_sdt <- "CR"
                  sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1
                }else{
                  iResponse_sdt <- "Hit"
                  sdt.mat_t[[L_x]][t,iResponse_sdt] <- sdt.mat_t[[L_x]][t,iResponse_sdt]+1}
              }
            }
            
            
          }
    
          providedResponses_recoded_t[[L_x]] <- c(unlist(providedResponses_recoded_t[[L_x]]),iResponse_recoded) 
        
        }
        
        
      }
      return(list(Scoring_Type1 = providedResponses_recoded_t, Scoring_Type2 = sdt.mat_t))
    }
    t_providedResponses_recoded <- recodeFx(correctLabels_per_level=correctLabels_per_level,
                                            w=t_matchesNomatches, x=t_labels,y=t_presentedStimuli,z=t_providedResponses)
    for(addToList in 1:3){
      providedResponses_recoded[[addToList]] <- rbind(providedResponses_recoded[[addToList]],
                                                      t_providedResponses_recoded$Scoring_Type1[[addToList]])
    }
    for(addToList in 1:3){
      sdt.mat[[addToList]] <- sdt.mat[[addToList]] + t_providedResponses_recoded$Scoring_Type2[[addToList]]
    }
    
  }
  
  for(addNameToList in 1:3){
    rownames(providedResponses_recoded[[addNameToList]]) <- paste("Judge_phase_",c(1:nJudgePhases),sep='')
    colnames(providedResponses_recoded[[addNameToList]]) <- paste("Q",c(1:nQuestions),sep='')
  }
  
  M_JS.phaseAndLevel <- lapply(providedResponses_recoded,function(i){rowMeans(matrix(unlist(i),ncol=6))-3})
  
  P_c_phaseAndLevel1 <- lapply(sdt.mat, function(i){
    f_correct <- rowSums(i[,c("Hit","CR")])
    f_wrong <- rowSums(i[,c("FA","Miss")])
    P_c <- f_correct/(f_correct+f_wrong)
    return(P_c)}
  )
  P_c_phaseAndLevel2 <- lapply(sdt.mat, function(i){
    f_correct <- rowSums(i[,c("Hit","CR")])
    P_c <- f_correct/rowSums(i)
    return(P_c)}
  )
  
  d_phaseAndLevel <- lapply(sdt.mat, function(i){
    i_relFreq <- prop.table(i,1)
    d_prime <- qnorm(pmax(0.001, pmin(i_relFreq[,"Hit"], 0.999)))-
      qnorm(pmax(0.001, pmin(i_relFreq[,"FA"], 0.999)))
    return(d_prime)
  })
  
  c_phaseAndLevel <- lapply(sdt.mat, function(i){
    i_relFreq <- prop.table(i,1)
    c_score <- -0.5*(pmax(0.001, pmin(i_relFreq[,"Hit"], 0.999)))+
      qnorm(pmax(0.001, pmin(i_relFreq[,"FA"], 0.999)))
    return(c_score)
  })
  
  d_df <- data.frame(d_prime=unlist(d_phaseAndLevel))
  c_score_df <- data.frame(c_score=unlist(c_phaseAndLevel))
  P_c1_df <- data.frame(P_c1=unlist(P_c_phaseAndLevel1))
  P_c2_df <- data.frame(P_c2=unlist(P_c_phaseAndLevel2))
  M_JS_df <- data.frame(M_JS=unlist(M_JS.phaseAndLevel))
  phases_and_levels <- strsplit(rownames(d_df),"", " ")
  phase <- unlist(lapply(phases_and_levels,function(i){paste("P",i[3],sep="")}))
  level <- unlist(lapply(phases_and_levels,function(i){paste("L",i[2],sep="")}))
  
  code.df <- data.frame(code, oSpan.score, WMC.group, taxo.type, phase,level,
                        d_prime=d_df,c_score=c_score_df,
                        P_c1=P_c1_df, P_c2=P_c2_df, M_JS=M_JS_df)
  
  return(code.df)
  
}}

fx_dinoNimi_judgment <- function(file, nRatingSteps, oSpan_Scores_df){
  
	data <- read.table(file,sep=",",fill=T,head = T)
	
	code <- c(as.matrix(data[3,"responses"]))
	code <- unlist(strsplit(code,":"))[2]
	code <- toupper(gsub("\\}|\\}|\"","",code))
	
	oSpan.score <- oSpan_Scores_df$nCorrect[oSpan_Scores_df$code==code]
	oSpan.group <- oSpan_Scores_df$WMC.group[oSpan_Scores_df$code==code]
	
	data_judgement <- data[which(data$summary_judgment =="summary_judgment"),]
	
	# Determining which exemplars were learned on category vs. item level
	presentedLabels <- as.matrix(data_judgement$presentedNamesLabels)[1]
	presentedLabels <- unlist(strsplit(presentedLabels,","))
	if(is.element("Edmontosaurus", presentedLabels)){taxo.type="taxo.dino"}else{taxo.type="taxo.ptero"}
	  
	timeStamps <- data_judgement[,"time_elapsed"]
	matchesNomatches <- as.matrix(data_judgement[,"matchNomatch"])
	matchesNomatches <- t(apply(matchesNomatches,1,function(i){as.numeric(unlist(strsplit(i,",")))}))
	providedResponses <- as.matrix(data_judgement[,"judgments"])
	providedResponses <- t(apply(providedResponses,1,function(i){as.numeric(unlist(strsplit(i,",")))}))
	
	presentedStimuli <- strsplit(as.matrix(data_judgement$presentedPicturesLabels),",")
	
	nJudgePhases <- nrow(providedResponses)
	nQuestions <- ncol(providedResponses)
	nRatingSteps <- 5
	# for loop to recode responses: whether the responses of e.g. "sure yes" or "sure no" are correct depends on
	# whether the presented picture and label have matched. A correct response provided with high confidence (e.g.
	# "sure yes" when picture and label have matched, is recoded as 5, a correct response with low confidence
	# as 4, and so on)
	
	# js: 0=sure yes, 4 = sure no
	
	for(t in 1: nJudgePhases){  # each row corresponds to one judgment phase
		if(t==1){providedResponses_recoded <- NULL}
		t_matchesNomatches <- matchesNomatches[t,]
		t_providedResponses <- providedResponses[t,]
		recodeFx <- function(x, y){
			for(i in 1: nQuestions){
				if(i==1){providedResponses_recoded <- NULL}
				iResponse <- y[i]
				iMatch <- x[i]
				if(iMatch==1){
				  recodeVec <- nRatingSteps:1
				  
				}else{
				  recodeVec <- 1:nRatingSteps
				  
				    
				    }
				iResponse_recoded <- recodeVec[iResponse+1]
				providedResponses_recoded <- c(providedResponses_recoded, iResponse_recoded)
			}
			return(providedResponses_recoded)
		}
		t_providedResponses_recoded <- recodeFx(t_matchesNomatches, t_providedResponses)
		providedResponses_recoded <- rbind(providedResponses_recoded, t_providedResponses_recoded)
	}
	"if(nJudgePhases!=8){
		nNAs <- rep(NA,(8-nJudgePhases)* nQuestions)
		NA_rows <- matrix(nNAs,nrow=(8-nJudgePhases))
		providedResponses_recoded <- rbind(providedResponses_recoded, NA_rows)
	}"
	rownames(providedResponses_recoded) <- paste("Judge_phase_",c(1:nrow(providedResponses_recoded)),sep='')
	colnames(providedResponses_recoded) <- paste("Q",c(1:nQuestions),sep='')
	
	"PhaseByResponse_FreqTab <- t(apply(providedResponses_recoded,1,function(x){
		probVec <- rep(0, nRatingSteps)
		probVec [as.numeric(names(table(x)))] <- table(x)
		return(probVec)
	}))
	colnames(PhaseByResponse_FreqTab) <- paste('P(',c(-2,-1,0,1,2),')',sep='')	
	if(nJudgePhases!=8){
		PhaseByResponse_FreqTab[(9-(8-nJudgePhases)):8,] <- NA
	}"
	
	average_response_perPhase <- rowMeans(providedResponses_recoded)
	
	presentedStimuli_sequence <- unlist(presentedStimuli)
	recodedResponses_sequence <- c(t(providedResponses_recoded))
	Level_of_learning <- rep("item_learning",length(presentedStimuli_sequence))
	if(taxo.type=="taxo.dino"){
		cat.episodes <- which(match(presentedStimuli_sequence,c("Pteranodontia","Rhamphorhynchinae"))>0)
		if(length(cat.episodes)>0){
			Level_of_learning[cat.episodes] <- "category_learning"
		}
		
	}else{
		cat.episodes <- which(match(presentedStimuli_sequence,c("Ornithoscelida","Saurolophidae"))>0)
		if(length(cat.episodes)>0){
			Level_of_learning[cat.episodes] <- "category_learning"
		}
	}
	phases_sequence <- c(sapply(c(1:length(average_response_perPhase)),function(i){rep(paste("P",i,sep=""),6)}))
	matchNomatch_seq <- c(t(matchesNomatches))
	"if(nJudgePhases!=8){
		Level_of_learning <- c(Level_of_learning, rep(NA,(8-nJudgePhases)*6))
		presentedStimuli_sequence <- c(presentedStimuli_sequence,rep(NA,(8-nJudgePhases)*6))
		matchNomatch_seq <- c(matchNomatch_seq,rep(NA,(8-nJudgePhases)*6))
	}"
	#ResponseByPhaseAndPerson <- data.frame(code, oSpan.score, WMC.group, taxo.type, phases_sequence, presentedStimuli_sequence, Level_of_learning, recodedResponses_sequence)
	#colnames(ResponseByPhaseAndPerson) <- c("code","oSpan","oSpan.group","taxo.type","phase","stimulus","LoL","response")
	ResponseByPhaseAndPerson <- data.frame(code, oSpan.score, oSpan.group, taxo.type, phases_sequence, presentedStimuli_sequence, Level_of_learning, recodedResponses_sequence)
	colnames(ResponseByPhaseAndPerson) <- c("code","oSpan","oSpan.group","taxo.type","phase","stimulus","LoL","response")
	
	# Hit and FA
	ResponseByPhaseAndPerson <- cbind(ResponseByPhaseAndPerson, matchNomatch_seq)
	l_Hits <- matchNomatch_seq==1&ResponseByPhaseAndPerson$response>=4
	l_Misses <- matchNomatch_seq==1&ResponseByPhaseAndPerson$response<=2
	l_CRs <- matchNomatch_seq==0&ResponseByPhaseAndPerson$response>=4
	l_FAs <- matchNomatch_seq==0&ResponseByPhaseAndPerson$response<=2
	l_dunno <- ResponseByPhaseAndPerson$response==3
	Hit.FA.vec <- rep(NA,nrow(ResponseByPhaseAndPerson))
	Hit.FA.vec[l_Hits] <- "Hit"
	Hit.FA.vec[l_Misses] <- "Miss"
	Hit.FA.vec[l_CRs] <- "CR"
	Hit.FA.vec[l_FAs] <- "FA"
	Hit.FA.vec[l_dunno] <- "Dunno"
	ResponseByPhaseAndPerson <- cbind(ResponseByPhaseAndPerson, Hit.FA.vec)
	Hit.vec <- rep(0,nrow(ResponseByPhaseAndPerson))
	Hit.vec[l_Hits] <- 1
	FA.vec <- rep(0,nrow(ResponseByPhaseAndPerson))
	FA.vec[l_FAs] <- 1
	CR.vec <- rep(0,nrow(ResponseByPhaseAndPerson))
	CR.vec[l_CRs] <- 1
	Miss.vec <- rep(0,nrow(ResponseByPhaseAndPerson))
	Miss.vec[l_Misses] <- 1
	Dunno.vec <- rep(0,nrow(ResponseByPhaseAndPerson))
	Dunno.vec[l_dunno] <- 1
	ResponseByPhaseAndPerson <- cbind(ResponseByPhaseAndPerson,Hit.vec,FA.vec,CR.vec,Miss.vec,Dunno.vec)
	
	#colnames(ResponseByPhaseAndPerson) <- c("code","oSpan","oSpan.group","taxo.type","phase","stimulus","LoL","response","matches","Hit.FA")
	colnames(ResponseByPhaseAndPerson) <- c("code","oSpan","oSpan.group","taxo.type","phase","stimulus","LoL","response","matches","Hit.FA","Hit","FA","CR","Miss","Dunno")
	
	toBeReturned <- list(code,ResponseByPhaseAndPerson)
	names(toBeReturned) <- c("code","ResponseByPhaseAndPerson")
	
	return(toBeReturned)
	
}

addUpMatrices <- function(x,y){
		x2 <- c(x)
		y2 <- c(y)
		addUp_i <- function(i,x2,y2){
			x2_i <- x2[i]
			y2_i <- y2[i]
			
			if(is.na(x2_i)&is.na(y2_i)){i_new <- NA}else if(is.na(x2_i)|is.na(y2_i)){
				i_new <- c(x2_i, y2_i)[!is.na(c(x2_i, y2_i))]}else{
					i_new <- x2_i+ y2_i}
			return(i_new)
			}
		addedUp <- sapply(c(1:length(x2)),x2=x2,y2=y2,addUp_i)		
		addedUp_matrix <- matrix(addedUp,nrow=8)	
		return(addedUp_matrix)
	}

ave.resp.fx <- function(low, splitCrit1, splitCrit2, phase, data){
	
	if(low==TRUE){rows.group <- data $WMC.group<= splitCrit1}else{
		rows.group <- data $WMC.group>= splitCrit2
	}
	
	rows.phase <- data $PhaseOfJudgment==phase
	group.by.phase <- data[ rows.group & rows.phase , ]
	
	group.by.phase.scores <- with(group.by.phase,c(mean(averageResponse,na.rm=T),sd(averageResponse,na.rm=T)))
	names(group.by.phase.scores) <- c("M","SD")
	return(group.by.phase.scores)
	
}

fx_dinoNimi_selection <- function(file, oSpan_Scores_df){
	
	#data <- read.table("118.CSV",sep=",",fill=T,head = T)
	#data <- read.table("429.CSV",sep=",",fill=T,head = T)
	#data <- read.table("418.CSV",sep=",",fill=T,head = T)
	data <- read.table(file,sep=",",fill=T,head = T)
		
	code <- c(as.matrix(data[3,"responses"]))
	code <- unlist(strsplit(code,":"))[2]
	code <- toupper(gsub("\\}|\\}|\"","",code))
	df_file_x <- data.frame(code)
	
	# Variable WMC
	oSpan.score <- oSpan_Scores_df$nCorrect[oSpan_Scores_df$code==code]
	WMC.group <- oSpan_Scores_df$WMC.group[oSpan_Scores_df$code==code]
	df_file_x <- cbind(df_file_x, oSpan.score, WMC.group)
	
	# Determining which exemplars were learned on category vs. item level
	displayed.buttons <- as.matrix(data[which(data$summary_judgment =="summary_judgment"),"presentedNamesLabels"])[1]
	displayed.buttons <- unlist(strsplit(displayed.buttons,","))
	if(is.element("Edmontosaurus", displayed.buttons)){taxo.type="taxo.dino"}else{taxo.type="taxo.ptero"}
	df_file_x <- cbind(df_file_x, taxo.type)
	
	data_categorization <- data[data$summary_categorization=="summary_categorization",]
	
	# Variable Time
	timeStamps <- data_categorization[,"time_elapsed"]
	viewingTimes <- data_categorization$rt
	catPhases <- data_categorization[,"genPhase"]
	df_file_x <- cbind(df_file_x, timeStamps, viewingTimes, catPhases)
	
	# Self-selected Stimuli
	displayedPics <- as.matrix(data_categorization[,"stimulus"])
	Stimuli <- c(as.matrix(sapply(displayedPics,function(x){gsub("img/|\\d|.png","",x)})))
	if(length(Stimuli =="Dorygnathusathus")>0){
		Stimuli[match("Dorygnathusathus", Stimuli)] <- "Dorygnathus"
	}
	N_SearchPaths <- sapply(c(1:length(Stimuli)), function(i,Stimuli){length(unique(Stimuli[1:i]))},Stimuli=Stimuli)
	requiredButtons <- c(as.matrix(data_categorization$requiredButton))
	N_Categories <- sapply(c(1:length(requiredButtons)), function(i, requiredButtons){length(unique(requiredButtons[1:i]))}, requiredButtons = requiredButtons)
	
	df_file_x <- cbind(df_file_x,Stimuli, N_SearchPaths, requiredButtons, N_Categories)
	
	# Variable Level of learning
	Level_of_learning <- rep("item_learning",length(Stimuli))
	if(taxo.type=="taxo.dino"){
		cat.episodes <- which(match(Stimuli,c("Pteranodon","Istiodactylus","Dorygnathus","Scaphognathus"))>0)
		if(length(cat.episodes)>0){
			Level_of_learning[cat.episodes] <- "category_learning"
		}
		
	}else{
		cat.episodes <- which(match(Stimuli,c("Parasaurolophus","Edmontosaurus","Rinchenia","Pisanosaurus"))>0)
		if(length(cat.episodes)>0){
			Level_of_learning[cat.episodes] <- "category_learning"
		}
	}
	df_file_x <- cbind(df_file_x, Level_of_learning)
	
	# Variables Practice & Exploration: ingredients
	
	nSelected_CurCat <- rep(NA,nrow(df_file_x))
	practiceAndExploreFx <- function(i, requiredButtons){
		CurCat <- requiredButtons[i]
		if(i==1){
			N <- 1
			dudi <- 0
			chi2 <- 0
			return(c(CurCat, 0, 0, 0, 0, 0, N,dudi, chi2))
			
		}else{
			categoriesUntilNow <- requiredButtons[1:i]
			freqs <- table(categoriesUntilNow)
			ps <- table(categoriesUntilNow)/sum(table(categoriesUntilNow))
			N <- length(ps)
			dudi <- sqrt(sum((ps-mean(ps))^2)/N)
			
			freqs_exp <- sum(freqs)/N
			chi2 <- sum((freqs-freqs_exp)^2/freqs_exp)
			
			previous.categories <- requiredButtons[1:(i-1)]
			previous.categories.tab <- table(previous.categories)
				
			if(is.element(CurCat, previous.categories)){
				
				#positions.curCat <- which(requiredButtons[1:i] ==CurCat)
				
				# NII = Number Intervening Items
				#NII.Fx <- function(pos.x, positions.curCat){ 
				#	if(pos.x==1){
				#			NII.x <- positions.curCat[pos.x]-1
				#		}else{
				#			NII.x <- (positions.curCat[pos.x]-positions.curCat[pos.x-1])-1
				#		}
				#	return(NII.x)
				#}
				#NII.vec <- sapply(c(1:length(positions.curCat)),NII.Fx, positions.curCat = positions.curCat)
				
				nCurCat <- previous.categories.tab[which(names(previous.categories.tab)==CurCat)]
				indices.CompCats <- which(names(previous.categories.tab)!=CurCat)
				nCompCats_unique <- length(indices.CompCats)
				if(length(indices.CompCats)>0){
					nCompCat <- sum(previous.categories.tab[indices.CompCats])
					#previous.p <- prop.table(previous.categories.tab[indices.CompCats])
					previous.p <- prop.table(previous.categories.tab)
					H <- -sum(previous.p *log2(previous.p))
					
					#denom <- log2(length(indices.CompCats))
					denom <- log2(length(previous.p))
					if(denom==0){
						H_norm <- 0
						practice_index <- 0
					}else{
						H_norm <- H/denom
						practice_index <- (nCurCat* (2-H_norm))/(nCurCat* (2-H_norm) +nCompCat*H_norm)
					}
				}else{
					nCompCat <- 0
					H_norm <- 0
					practice_index <- (nCurCat* (2-H_norm))/(nCurCat* (2-H_norm) +nCompCat*H_norm)
				}
			}else{
				nCurCat <- 0
				indices.CompCats <- which(names(previous.categories.tab)!=CurCat)
				nCompCats_unique <- length(indices.CompCats)
				if(length(indices.CompCats)>0){
					nCompCat <- sum(previous.categories.tab[indices.CompCats])
					#previous.p <- prop.table(previous.categories.tab[indices.CompCats])
					previous.p <- prop.table(previous.categories.tab)
					H <- -sum(previous.p *log2(previous.p))
					#denom <- log2(length(indices.CompCats))
					denom <- log2(length(previous.p))
					if(denom==0){
						H_norm <- 0
						practice_index <- 0
					}else{
						H_norm <- H/denom
						practice_index <- 0
					}
				}else{
					nCompCat <- 0
					H_norm <- 0
					practice_index <- 0
				}
			}
			
			return(c(CurCat, nCurCat, nCompCat, nCompCats_unique, H_norm, practice_index,N,dudi, chi2))
			
		}
		
	}
	practiceAndExplore.df <- data.frame(t(sapply(c(1:length(requiredButtons)),practiceAndExploreFx, requiredButtons=requiredButtons)))
	colnames(practiceAndExplore.df) <- c("CurCat", "nCurCat", "nCompCat", "nCompCats_unique", "H_norm", "practice_index","N","dudi","chisq")
	
	df_file_x <- cbind(df_file_x, practiceAndExplore.df[,c("nCurCat", "nCompCat", "nCompCats_unique", "H_norm", "practice_index","N","dudi","chisq")])
	
	chosenButtons <- data_categorization$chosenButton
	matchNomatch <- data_categorization$correctCategorization
		
	for(i in 2:nrow(df_file_x)){
	  if(i==2){rep.df <- data.frame(rep.i=0,lag.i=NA)}
	  
	  if(is.element(df_file_x$Stimuli[i],df_file_x$Stimuli[1:(i-1)])){
	    rep.i <- 1
	    lag.i <- i-max(which(df_file_x$Stimuli[i]==df_file_x$Stimuli[1:(i-1)]))
	  }else{
	    rep.i <- 0
	    lag.i <- NA
	  }
	  rep.df <- rbind(rep.df,data.frame(rep.i,lag.i))
	}
	
	df_file_x <- cbind(df_file_x, chosenButtons, matchNomatch, rep.df, i=1:nrow(df_file_x) )
	
	stimuli.code.x <- c(as.matrix(df_file_x $Stimuli))
	switchFx <- function(i, stimuli.code.x){
		if(length(unique(c(stimuli.code.x[i],stimuli.code.x[i-1])))==2){
			switch.i <- "switch"
		}else{
			switch.i <-"noSwitch"
		}
		return(switch.i)
	}
	switch.vec <- c(NA,sapply(c(2:length(stimuli.code.x)),switchFx, stimuli.code.x= stimuli.code.x))
	df_file_x <- cbind(df_file_x,switch.vec)
	
	## Binary stimulus coding
	r.sumCat <- which(data$summary_categorization=="summary_categorization")
	f1 <- data[r.sumCat-4,"button_pressed"]
	f2 <- data[r.sumCat-3,"button_pressed"]
	f3 <- data[r.sumCat-2,"button_pressed"]
	name.subordinate <- data[r.sumCat-1,"stimulus"]
	name.subordinate <- c(as.matrix(sapply(name.subordinate,function(x){gsub("img/|\\d|.png","",x)})))
	if(length(name.subordinate =="Dorygnathusathus")>0){
	  name.subordinate[match("Dorygnathusathus", name.subordinate)] <- "Dorygnathus"
	}
	taxo.subordinate <- c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus",
	                      "Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
	names(taxo.subordinate) <- c("Ornithoscelida","Ornithoscelida","Saurolophidae","Saurolophidae",
	                             "Rhamphorhynchinae","Rhamphorhynchinae","Pteranodontia","Pteranodontia")
	name.basicLevel <- names(taxo.subordinate)[match(name.subordinate,taxo.subordinate)]
	df_file_x <- cbind(df_file_x,f1,f2,f3,name.basicLevel,name.subordinate)
	
	return(df_file_x)
	
}


## get data for simulation
fx_featurePatterns <- function(file, oSpan_Scores_df){
  
  #data <- read.table("118.CSV",sep=",",fill=T,head = T)
  #data <- read.table("429.CSV",sep=",",fill=T,head = T)
  #data <- read.table("418.CSV",sep=",",fill=T,head = T)
  data <- read.table(file,sep=",",fill=T,head = T)
  
  code <- c(as.matrix(data[3,"responses"]))
  code <- unlist(strsplit(code,":"))[2]
  code <- toupper(gsub("\\}|\\}|\"","",code))
  df_file_x <- data.frame(code)
  
  # Variable WMC
  oSpan.score <- oSpan_Scores_df$nCorrect[oSpan_Scores_df$code==code]
  WMC.group <- oSpan_Scores_df$WMC.group[oSpan_Scores_df$code==code]
  df_file_x <- cbind(code, oSpan.score, WMC.group)
  
  # Determining which exemplars were learned on category vs. item level
  displayed.buttons <- as.matrix(data[which(data$summary_judgment =="summary_judgment"),"presentedNamesLabels"])[1]
  displayed.buttons <- unlist(strsplit(displayed.buttons,","))
  if(is.element("Edmontosaurus", displayed.buttons)){taxo.type="taxo.dino"}else{taxo.type="taxo.ptero"}
  df_file_x <- cbind(df_file_x, taxo.type)
  
  r.sumCat <- which(data$summary_categorization=="summary_categorization")
  f1 <- data[r.sumCat-4,"button_pressed"]
  f2 <- data[r.sumCat-3,"button_pressed"]
  f3 <- data[r.sumCat-2,"button_pressed"]
  name.subordinate <- data[r.sumCat-1,"stimulus"]
  name.subordinate <- c(as.matrix(sapply(name.subordinate,function(x){gsub("img/|\\d|.png","",x)})))
  if(length(name.subordinate =="Dorygnathusathus")>0){
    name.subordinate[match("Dorygnathusathus", name.subordinate)] <- "Dorygnathus"
  }
  requiredButton <- c(as.matrix(data[r.sumCat,"requiredButton"]))
  chosenButton <- c(as.matrix(data[r.sumCat,"chosenButton"]))
  rightOrWrong <- data[r.sumCat,"correctCategorization"]
  
  taxo.subordinate <- c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus",
                                "Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")
  names(taxo.subordinate) <- c("Ornithoscelida","Ornithoscelida","Saurolophidae","Saurolophidae",
                                       "Rhamphorhynchinae","Rhamphorhynchinae","Pteranodontia","Pteranodontia")
  
  name.basicLevel <- names(taxo.subordinate)[match(name.subordinate,taxo.subordinate)]
  searchPhase <- data[r.sumCat,"genPhase"]
  df_file_x <- cbind(df_file_x,
          data.frame(searchPhase,f1,f2,f3,name.subordinate,name.basicLevel,requiredButton,chosenButton,rightOrWrong)
        )
  
  return(df_file_x)
  
}


