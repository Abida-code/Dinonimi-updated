##### modeling dinoNimi
rm(list=ls())
setwd("~/Dropbox/Projekte/CEITER/Sea&Sense/SUSTAIN_Modeling")

## VARIABLES' meaning
# I ... stimulus vector
# pos ... position of stimulus in the multidimensional representational space
# i ... ith stimulus dimension
# k ... kth feature on dimension i
# j ... cluster index 
# H_j_act ... activation of jth cluster H
# H_j_out ... output of jth cluster after lateral competition 
# m ... number of stimulus dimensions
# V ... vector representing each number vi of features per stimulus dimension i
# lambda_i ... tuning of the receptive field for the ith input dimension
# z ... queried stimulus dimension
# C_out_zk ... output of the output unit representing the kth value of z

## PARAMETERS' meaning
# r ... attentional parameter (accentuates the effect of the lambda vector)
# beta ... cluster competition
# eta ... learning rate
# d ... response parameter ("When d is high, accuracy is stressed and the output unit with the largest output is almost always chosen", p. 315)

## FUNCTIONS
source("SUSTAIN_fx.R")

setwd("~/Dropbox/Projekte/CEITER/Sea&Sense/SUSTAIN_Modeling")
source("prepare_emp.R")

## parameter setting

par.x <- c(2.844642,2.386305,12.0,0.09361126)
names(par.x) <- c("r","beta","d","eta")

## stimulus properties
labelQuery <- TRUE
m <- 4
V <- c(2,2,2,6)
dino.structure <- list(c(1,1,1,1),
                       c(1,1,2,2),
                       c(1,2,1,3),
                       c(1,2,2,4),
                       c(2,1,1,5),
                       c(2,1,2,5),
                       c(2,2,1,6),
                       c(2,2,2,6))
ptero.structure <- list(c(1,1,1,1),
                        c(1,1,2,1),
                        c(1,2,1,2),
                        c(1,2,2,2),
                        c(2,1,1,3),
                        c(2,1,2,4),
                        c(2,2,1,5),
                        c(2,2,2,6))
#stimulus_structure <- shepard_structure_I
structure_name <- "Type1"
feature_coding <- c(0,1)

## Srun settings
ResponseSet.dinoTaxo <- c("Pisanosaurus","Rinchenia","Edmontosaurus","Parasaurolophus",
                          "Rhamphorhynchinae","Pteranodontia")
ResponseSet.pteroTaxo <- c("Ornithoscelida","Saurolophidae",
                           "Dorygnathus","Scaphognathus","Istiodactylus","Pteranodon")


y <- 0
nRun <- 2
repeat{
  y <- y+1
  if(y==1){
    performance.df <- NULL
    mentalStructure.df <- NULL
  }
  for(x in 1:length(listedFilesDinoNimi)){
    if(x==1){t <- Sys.time()}
    
    file.x <- listedFilesDinoNimi[x]
    data_student <- fx_student_data(file.x, oSpan_Scores_df)
    code <- c(as.matrix(data_student$code[1]))
    group.x <- as.numeric(c(as.matrix(data_student$WMC.group[1])))
    if((group.x==1)|(group.x==4)){
      nTrials <- nrow(data_student)
      structure_name <- c(as.matrix(data_student$taxo.type[1]))
      selfSelectedStimuli <- data_student[,c("f1","f2","f3")]+1
      selfSelectedStimuli_strings <- apply(selfSelectedStimuli,1,function(i){toString(i)})
      
      if(structure_name=="taxo.dino"){
        stimulus_structure <- dino.structure  
        ResponseSet <- ResponseSet.dinoTaxo
      }else{
        stimulus_structure <- ptero.structure
        ResponseSet <- ResponseSet.pteroTaxo
      }
      for(z in 1:8){
        if(z==1){selfSelectedStimuli_trueCats <- rep(NA,nTrials)}
        j <- toString(stimulus_structure[[z]][1:3])
        r_j <- which(selfSelectedStimuli_strings==j)
        if(length(r_j)>0){selfSelectedStimuli_trueCats[r_j] <- stimulus_structure[[z]][4]}
      }
      
      selfSelectedStimuli_chosenCats <- match(data_student$chosenButton,ResponseSet)
      selfSelectedStimuli <- cbind(data_student$searchPhase,selfSelectedStimuli,selfSelectedStimuli_trueCats,selfSelectedStimuli_chosenCats,data_student$rightOrWrong)
      
      colnames(selfSelectedStimuli) <- c("searchPhase","f_sup","f_basic","f_sub","f_cat_true","f_cat_obs","RightOrWrong")
      
      agent.outcome <- agent_fx(par.x=par.x, data_student=data_student,
                                labelQuery, structure_name, stimulus_structure, selfSelectedStimuli, ResponseSet,
                                m, V)
      
      performance.df <- rbind(performance.df,data.frame(y,code,group.x,agent.outcome$session_df))
      
      agent.lambdas <- agent.outcome$lambda_vec
      mentalStructure.agent <- data.frame(y,code,group.x,agent.lambdas[1],agent.lambdas[2],agent.lambdas[3],nrow(agent.outcome$H))
      names(mentalStructure.agent) <- c("y","code","l_sup","l_basic","l_sub","#clusters")
      mentalStructure.df <- rbind(mentalStructure.df,mentalStructure.agent) 
    }
    print(c(y,paste(x," out of ",length(listedFilesDinoNimi) ),
            paste("Seconds: ",round(c(as.matrix(Sys.time()-t)),2))))
  }
  y.duration <- c(as.matrix(Sys.time()-t))
  if(y==1){duration.vec <- y.duration}else{duration.vec <- c(duration.vec,y.duration)}
  if(y==nRun){break}
  
}
mean(duration.vec)

mean(mentalStructure.df[,"#clusters"])
colMeans(mentalStructure.df[,3:5])

P_E_observed_overall <- aggregate((1-correctResponse_student) ~ trial, performance.df, function(i){c(mean(i,na.rm=T),sd(i,na.rm=T))})
P_E_predicted_overall <- aggregate((1-correctResponse_sustain) ~ trial, performance.df, function(i){c(mean(i,na.rm=T))})
P_E_observed <- aggregate((1-correctResponse_student) ~ trial*group.x, performance.df, function(i){c(mean(i,na.rm=T),sd(i,na.rm=T))})
P_E_predicted <- aggregate((1-correctResponse_sustain) ~ trial*group.x, performance.df, function(i){c(mean(i,na.rm=T))})

nReps <- 17 #30
obs.M <- P_E_observed_overall[,2][1:nReps,1]
obs.SD <- P_E_observed_overall[,2][1:nReps,2]
n <- length(unique(performance.df$code))
Error.lO.Conf95 <- qnorm(0.975)* obs.SD/sqrt(n)

plot(obs.M,ylim=c(0,1),type="p",lwd=.5,pch=16,xlab=expression("Practice"*italic(" n")),ylab=expression(italic("P(E)")*" = Prob. choosing incorrect label"),cex.axis=.8)
points(P_E_predicted_overall[1:nReps,2],col="green",pch=16)

trials <- 1:nReps
points(P_E_predicted_overall[1:nReps,2],pch=1,col="blue")

sustain.M <- P_E_predicted_overall[1:nReps,2]
nls.fit.sustain <- nls(y~a+(1-x^-b),data=data.frame(x=c(1:length(sustain.M)),y= sustain.M),start=list(a=.2,b=.5))
interc.sustain <- summary(nls.fit.sustain)$coefficients[1]
decay.sustain <- summary(nls.fit.sustain)$coefficients[2]
nls.pred.sustain <- interc.sustain + (1-c(1:nReps)^-decay.sustain)
nls.pred.toPlot.sustain <- interc.sustain + (1-seq(1, nReps,length=100)^-decay.sustain)
#points(nls.pred.toPlot.sustain ~seq(1, nReps,length=100),type="l",lwd=.5)
summary(lm(nls.pred.sustain ~obs.M))[9]

nls.fit.overall <- nls(y~a+(1-x^-b),data=data.frame(x=c(1:length(obs.M)),y= obs.M),start=list(a=.2,b=.5))
interc.overall <- summary(nls.fit.overall)$coefficients[1]
decay.overall <- summary(nls.fit.overall)$coefficients[2]
model.pred.overall <- interc.overall + (1-c(1:nReps)^-decay.overall)
model.pred.toPlot.overall <- interc.overall + (1-seq(1, nReps,length=100)^-decay.overall)
points(model.pred.toPlot.overall ~seq(1, nReps,length=100),type="l",lwd=.5)
summary(lm(model.pred.overall ~obs.M))[9]

legend(20,0.9,c("Observed","SUSTAIN"),pch=c(16,1),box.lwd = 0)

unique.codes <- unique(performance.df$code)
for(cx in 1:length(unique.codes)){
  if(cx==1){NIE.df <- NULL}
  code <- unique.codes[cx]
  r.code <- which(performance.df$code==code)
  code.data <- performance.df[r.code,]
  code.cats <- sort(unique(code.data$cat_true))
  
  for(x in 1:length(code.cats)){
    cat.x <- code.cats[x]
    r.cat.x <- which(code.data$cat_true==cat.x)
    if(length(r.cat.x)>1){
      cat.x.NIE <- sapply(c(2:length(r.cat.x)),function(i){r.cat.x[i]-r.cat.x[i-1]})
      r.NIE <- r.cat.x[2:length(r.cat.x)]
      cat.x.df <- data.frame(code,cat.x,cat.x.NIE,code.data[r.NIE,]$correctResponse_student,code.data[r.NIE,]$correctResponse_sustain)
      names(cat.x.df) <- c("code","cat","NIE","Correct_student","Correct_sustain")
      
      NIE.df <- rbind(NIE.df,cat.x.df)
    }
  }
}

nObs <- 10
E.NIE.student <- aggregate((1-Correct_student)~NIE,NIE.df,function(i){c(mean(i),sd(i))})
E.NIE.student <- E.NIE.student[1:nObs,]
E.NIE.student.M <- E.NIE.student[,2][,1]
E.NIE.student.SD <- E.NIE.student[,2][,2]
n <- length(unique.codes)
E.NIE.conf <- qnorm(0.975)* E.NIE.student.SD /sqrt(n)
E.NIE.sustain <- aggregate((1-Correct_sustain)~NIE,NIE.df,mean)
E.NIE.sustain.M <- E.NIE.sustain[1:nObs,2]

par(mar=c(4,4,1,4))
plot(E.NIE.student.M, ylim=c(0,.8), xlim=c(1,nObs),type="o",lwd=.5,pch=16,xlab=expression(italic("NIE")*" = Number Intervening Episodes"),ylab=expression(italic("P(E)")*" = Prob. choosing incorrect label"))
xs <- 1:nObs
arrows(xs, E.NIE.student.M + E.NIE.conf, xs, E.NIE.student.M-E.NIE.conf,lwd=.5,code=3,angle=90,length=.05)
points(E.NIE.sustain.M,type="o",lty=2,pch=1,col="green")

legend(1,0.8,box.lwd=0,pch=c(16,1),lty=c(1,2),c("Observed","Sustain"))



