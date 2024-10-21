### setting up data frames to get empirical data patterns

# read in some functions 
setwd("~/Dropbox/Class8c_7May")
source("fx_manifestAna.R")

# read in codes of al participants with valid data; can be read in from either the oSpan_pruned- or dinoNimi_pruned-folder
setwd("~/Dropbox/Class8c_7May/oSpan_pruned")
listedFilesOSpan <- list.files(pattern="*.CSV|*.txt")
validCodes <- toupper(c(as.matrix(sapply(listedFilesOSpan, code_fx))))

## creating a data.frame with all oSpan scores and corresponding codes
setwd("~/Dropbox/Class8c_7May/oSpan_pruned")
for(i in 1:length(listedFilesOSpan)){
  if(i==1){oSpan_Scores_df <- NULL}
  iData <- fxOspan(listedFilesOSpan[i])
  oSpan_Scores_df <- rbind(oSpan_Scores_df, iData)
}

solve_accuracies <- oSpan_Scores_df$solve_accuracy
recall_accuracy <- oSpan_Scores_df$nCorrect

# excluding outliers, i.e., participants whose solve-accuracy is below the average accuracy minus two times the SD
outliers <- which ( solve_accuracies <= mean(solve_accuracies) - 2*sd(solve_accuracies) )

WMC.group <- with(oSpan_Scores_df, cut(nCorrect, breaks=quantile(nCorrect), include.lowest=TRUE))
WMC.group <- as.numeric(WMC.group)
oSpan_Scores_df <- cbind(oSpan_Scores_df, WMC.group)

accu.outlier <- rep(TRUE,nrow(oSpan_Scores_df))
accu.outlier[outliers] <- FALSE
oSpan_Scores_df <- cbind(oSpan_Scores_df, accu.outlier)

### analyzing dinoNimi
setwd("~/Dropbox/Class8c_7May/dinoNimi_pruned")
listedFilesDinoNimi <- list.files(pattern="*.CSV|*.txt")
codes_dinoNimi <- toupper(c(as.matrix(sapply(listedFilesDinoNimi, code_fx))))

## judgments

for(file_x in 1:length(listedFilesDinoNimi)){
  if(file_x==1){
    selection.df <- NULL
  }
  file_x_selectionResults <- fx_dinoNimi_selection(listedFilesDinoNimi[file_x], oSpan_Scores_df= oSpan_Scores_df)
  selection.df <- rbind(selection.df, file_x_selectionResults)		
}
selection.df [which(selection.df $Stimuli=="Dorygnathusathus"),"Stimuli"] <- "Dorygnathus"
selection.df $N_SearchPaths <- as.numeric(as.matrix(selection.df $N_SearchPaths))
selection.df $N_Categories <- as.numeric(as.matrix(selection.df $N_Categories))
selection.df $nCurCat <- as.numeric(as.matrix(selection.df $nCurCat))
selection.df $nCompCat <- as.numeric(as.matrix(selection.df $nCompCat))
selection.df $nCompCats_unique <- as.numeric(as.matrix(selection.df $nCompCats_unique))
selection.df $H_norm <- as.numeric(as.matrix(selection.df $H_norm))
selection.df $practice_index <- as.numeric(as.matrix(selection.df $practice_index))
selection.df $matchNomatch <- as.numeric(as.matrix(selection.df $matchNomatch))
selection.df $timeStamps <- round(selection.df$timeStamps/(1000*60))
selection.df $code <- as.factor(selection.df $code)
selection.df $Level_of_learning <- as.factor(selection.df $Level_of_learning)
selection.df $N <- as.factor(selection.df $N)
selection.df $dudi <- as.factor(selection.df $dudi)
selection.df $catPhases <- as.factor(selection.df $catPhases)

## get data for simulation
# file <- listedFilesDinoNimi[1]
fx_student_data <- function(file, oSpan_Scores_df){
  
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


