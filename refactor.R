#Set working directory
path = "/home/caitlin/LabProjects/"
#Here is where packages and ggplot set-up chunk should go
setwd(path)
library(dplyr)
#TO-DO: relative paths

#SIMULATED DATA

#Simulation Data Structures

#NEW TYPES FOR CMRLOG
setClass("TrackedBird", representation(ID = "numeric", yearFirstCaptured = "numeric",
                                       stateList = "list"))

setClass("CMRState", representation(time = "numeric", wasSeen = "logical",
                                    wasCaptured = "logical", currentSite = "numeric",
                                    isSeropositive = "logical")) 
  
#Creates a naive individual. Stored as vector.
newIndividual = function(animalNum, #Numeric -- ID number for a given animal
                         siteNum){ #Numeric -- site the animal is recruited to
  individual = list(animalNum, siteNum, FALSE, FALSE, FALSE)
  names(individual) = c("animalNumber",
                 "currentSite",
                 "isExposed",
                 "isSeropositive",
                 "isMarked")
  return(as.data.frame(individual))
}


#Creates a vector that stores current number of animals in each state at
#a given site.
newSiteState = function(sitePop, #Data frame -- row vecs are individuals at site
                        siteNum, #Numeric -- ID of site that is being recorded
                        year){ #Numeric -- timestep number.
  popSize = nrow(sitePop)
  numSeropositive = sum(sitePop$isSeropositive == 1)
  siteState = list(year, siteNum, popSize, numSeropositive)
  names(siteState) = c("year",
                       "siteNumber",
                       "populationSize",
                       "numSeropositive")
  return(as.data.frame(siteState))
}


#Creates naive population. Population is data frame with rows of type Individual.
newPopulation = function(popSizeVec){ #Length of vector is number of sites.
                                      #Vector entries are pop size at each site.
  population = c()
  popCounted = 0 #For indexing animals
  for(i in 1:length(popSizeVec)){
    sampleSitePop = addIndividuals(popSizeVec[i], popCounted + 1, i)
    population = rbind(population, sampleSitePop)
    popCounted = popCounted + popSizeVec[i]
  }
  return(as.data.frame(population))
}


#Initializes vector of parameters for the population for one site
newEpiSiteVec = function(siteNum, #Site number
                         popSize, #Population at site
                         deathChance, #Probability that an animal dies
                         seroRevChance, #Prob that non-dying seropos animal reverts
                         forceOfInf){ #Prob that non-dying animal seroconverts
  epiSiteVec = list(siteNum, popSize, deathChance, seroRevChance, forceOfInf)
  names(epiSiteVec) = c("siteNumber", "popSize", "deathChance", "seroRevChance",
                        "infectionChance")
  return(as.data.frame(epiSiteVec))
}


#Internal screaming
newMarkSiteVec = function(siteNum, #Site number
                          numMarked, #Number at site that are sampled each year
                          resightProb, #Prob of marked animal being resighted
                          captureProb){ #Prob of resighted animal being captured
  markSiteVec = list(siteNum, numMarked, resightProb, captureProb)
  names(markSiteVec) = c("siteNumber", "numberMarked", "resightingChance",
                         "recaptureChance")
  return(as.data.frame(markSiteVec))
}
 

#Global capture mark recapture data log

#Global because without references, I otherwise cannot modify a value without
#returning it, so I have to alter the scope. Index names will be animal numbers. 
#Contents will be lists of type TrackedBird
CMRLOG = list()

#Global most recent ID to keep track of the IDs we've used
newestAnim = 0

#Methods


#Makes new individuals to be added to site. Each individual is numbered,
#starting with passed startingIndex value.
addIndividuals = function(numIndividuals, #numeric
                          startingIndex, #numeric
                          siteNum) { #numeric
  lastIndex = startingIndex + numIndividuals - 1
  popArray = c()
  for(i in startingIndex:lastIndex){
    thisAnimal = newIndividual(i, siteNum)
    popArray = rbind(popArray, thisAnimal)
  }
  newestAnim <<- lastIndex
  return(popArray)
}


#counts the number of individuals in each state and returns that as an array
#one row per site
getPopStats = function(population,#numeric
                       numSites, #numeric
                       year){ #numeric
  popStats = c()
  for(i in 1:numSites){
    sitePop = population[which(population$currentSite == i),]
    state = newSiteState(sitePop, i, year)
    popStats = rbind(popStats, state)
  }
  return(popStats)
}


#Given a marked individual, returns the state history of that individual
getIndividualHistory = function(ind){ #numeric
  if(ind$isMarked == FALSE){
    stop("Tried to get CMR history of unmarked individual.")
  }
}


#Rounds proportion of individuals to a whole number
roundProp = function(proportion, #numeric between 0 and 1
                     number){ #numeric
  return(round(number * proportion))
}


#gets pop at a site -- mainly here cuz i gotta do it so much
getSitePop = function(siteNum, population){
  sitePop = population[population$currentSite == siteNum,]
  return(sitePop)
}

#In bird system, dispersers will die with chance from their pre-dispersal
#site but become infected as a function of the force of infection at their
#post-dispersal site. Move function calls as necessary for specific system.
dispersalTransitions = function(epiSiteVec, #type EpiSiteVec
                                population, #type Population
                                siteNum, #numeric
                                dispersalVec){ #atomic vector
  sitePop = population[population$currentSite == siteNum,]
  sitePop = removeDying(epiSiteVec, sitePop)
  sitePop = disperse(sitePop, dispersalVec, epiSiteVec)
  return(sitePop)
}


#Removes dying individuals from the population at a site given a vec of model
#params
removeDying = function(epiSiteVec, #type EpiSiteVec
                       sitePopulation){ #type Population
  numReplaced = rbinom(1, nrow(sitePopulation), epiSiteVec$deathChance)
  dyingIndividuals = sample(sitePopulation[,"animalNumber"], numReplaced)
  sitePopulation = sitePopulation[!is.element(sitePopulation$animalNumber,
                                              dyingIndividuals),]
  return(sitePopulation)
}


#Does dispersal for a site
disperse = function(sitePop, #type Population
                    dispersalVec, #atomic vector
                    epiSiteVec){ #type EpiSiteVec
  sitePop$currentSite =sample(x = 1:length(dispersalVec), size = nrow(sitePop),
                              replace = TRUE, prob = dispersalVec)
  return(sitePop)
}


#Does births and infections which both happen after dispersal
postDispersalTransitions = function(epiSiteVec, #type EpiSiteVec
                                    population, #type Population
                                    dispersalVec){ #atomic vector
  newestID = newestAnim
  sitePop = getSitePop(epiSiteVec$siteNumber, population)
  sitePop = popReplacement(epiSiteVec, sitePop, newestID)
  sitePop = seroTransitions(epiSiteVec, sitePop)
  return(sitePop)
}


#Does births
popReplacement = function(epiSiteVec, #type EpiSiteVec
                          sitePopulation, #type Population
                          newestAnimalID){ #numeric - highest index in total pop
  numReplaced = epiSiteVec$popSize - nrow(sitePopulation)
  if(numReplaced > 0){
    newBirths = addIndividuals(numReplaced, newestAnimalID + 1,
                             epiSiteVec$siteNumber)
    sitePopulation = rbind(sitePopulation, newBirths)
  }
  return(sitePopulation)
}


#Does seroreversion and seroconversion for one site
seroTransitions = function(epiSiteVec, #type EpiSiteVec
                           sitePop){ #type Population
  seroposList = sitePop[sitePop$isSeropositive == TRUE,]
  numSerorev = rbinom(1, nrow(seroposList), epiSiteVec$seroRevChance)
  serorevList = sample(seroposList[,"animalNumber"], numSerorev)
  susceptibleList = sitePop[sitePop$isSeropositive == FALSE,]
  numSeroconv = rbinom(1, nrow(susceptibleList), epiSiteVec$infectionChance)
  seroconvList = sample(sitePop[,"animalNumber"], numSeroconv)
  sitePop[is.element(sitePop$animalNumber, serorevList),"isSeropositive"] = FALSE
  sitePop[is.element(sitePop$animalNumber, seroconvList), "isSeropositive"] = TRUE
  return(sitePop)
}

#Does longitudinal sampling for a given population state
captureMarkRecapture = function(population, #type Population
                               cmrMatrix, #data frame of type MarkSiteVec
                               time){ #simulation year
  updatedPop = c()
  numSites = nrow(cmrMatrix)
  for (i in 1:numSites){
    cmrVec = cmrMatrix[i,]
    sitePop = getSitePop(cmrVec$siteNumber, population)
    sitePop = cmrAtSite(sitePop, cmrVec, time)
    updatedPop = rbind(updatedPop, sitePop)
  }
  updateLogDead(time)
  return(updatedPop)
}


#Handles resighting and recapturing of already marked animals at one site
cmrAtSite = function(sitePop, #type Population
                     cmrVec, #type MarkSiteVec
                     time){ #simulation year
  markedInds = sitePop[sitePop$isMarked == TRUE,]
  unmarkedInds = sitePop[sitePop$isMarked == FALSE,]
  numSeen = rbinom(1, length(markedInds$animalNumber), cmrVec$resightingChance)
  numRecaptured = rbinom(1, numSeen, cmrVec$recaptureChance)
  seeAndCapture(markedInds, numSeen, numRecaptured, cmrVec$siteNumber, time)
  newMarkedInds = markNewIndividuals(unmarkedInds, cmrVec$numberMarked - numRecaptured,
                                     cmrVec$siteNumber, time)
  if(length(newMarkedInds$animalNumber) != 0) {
    sitePop[is.element(sitePop$animalNumber, newMarkedInds$animalNumber),]$isMarked = TRUE
  }
  return(sitePop)
}


#Handles the animals that are already marked and updates the CMR data log
seeAndCapture = function(marked, numSeen, numRecaptured, site, time){
  if(numSeen > 0){
    seenInds = sample_n(marked, numSeen)
    recapturedInds = sample_n(seenInds, numRecaptured)
    seenOnlyInds = filter(seenInds, !seenInds$animalNumber %in% recapturedInds$animalNumber)
    unseenInds = filter(marked, !marked$animalNumber %in% seenInds$animalNumber)
    updateLogSeen(seenOnlyInds, site, time)
    updateLogRecap(recapturedInds, site, time)
    updateLogNotSeen(unseenInds, site, time)
  }
}


#Marks new individuals at a site to make up for any deficits
#Updates CMR data log
markNewIndividuals = function(unmarkedInds,
                              numNeeded, siteNum, time){
  if(ncol(unmarkedInds)!=5){
    print("wgat?")
  }
  if(numNeeded > length(unmarkedInds$animalNumber)){
    updateLogNew(unmarkedInds, siteNum, time)
    return(unmarkedInds)
  }
  if(numNeeded < 0){ #we've caught more than expected due to dispersal!
    dummyDF = data.frame(animalNumber = integer())
    #Returns an empty named dataframe so the program will see there are
    #no animals in returnValue$animalNumber.
    return(dummyDF)
  }
  newMarkeds = sample_n(unmarkedInds, numNeeded)
  updateLogNew(newMarkeds, siteNum, time)
  return(newMarkeds)
}


#Updates CMR log to add new individuals
updateLogNew = function(newMarkeds, siteNum, t){
  if(length(newMarkeds$animalNumber) > 0){
    birdList = list()
    for(i in 1:length(newMarkeds$animalNumber)){
      state = new("CMRState", time = t, wasSeen = TRUE, wasCaptured = TRUE,
       currentSite = siteNum, isSeropositive = newMarkeds$isSeropositive[i])
      cmrList = list(state)
      capturedBird = new("TrackedBird", ID = newMarkeds$animalNumber[i],
                         yearFirstCaptured = t, stateList = cmrList)
      birdList[i] = capturedBird
    }
    names(birdList) = newMarkeds$animalNumber
    CMRLOG <<- append(CMRLOG,birdList)
  }
}


#Updates CMR log to chart which individuals have been seen
updateLogSeen = function(seenList, siteNum, t){
  if(length(seenList$animalNumber) > 0){
    for(i in 1:length(seenList$animalNumber)){
      seenState = new("CMRState", time = t, wasSeen = TRUE, wasCaptured = FALSE,
                  currentSite = siteNum, isSeropositive = NA)
      id = as.character(seenList$animalNumber[i])
      change = append(CMRLOG[[id]]@stateList,seenState)
      CMRLOG[[id]]@stateList <<- change
    }
  }
}


#Updates CMR log to chart which individuals have been recaptured
updateLogRecap = function(recapList, siteNum, t){
  if (length(recapList$animalNumber) > 0){
    for(i in 1:length(recapList$animalNumber)){
      recapState = new("CMRState", time = t, wasSeen = TRUE, wasCaptured = TRUE,
                  currentSite = siteNum, isSeropositive = recapList$isSeropositive[i])
      id = as.character(recapList$animalNumber[i])
      change = append(CMRLOG[[id]]@stateList,recapState)
      CMRLOG[[id]]@stateList <<-change
    }
  }
}


#Updates CMR log to record unseen individuals
updateLogNotSeen = function(unseenList, siteNum, t){
  if(length(unseenList$animalNumber > 0)){
    for(i in 1:length(unseenList$animalNumber)){
      unseenState = new("CMRState", time = t, wasSeen = FALSE, wasCaptured = FALSE,
                  currentSite = 0, isSeropositive = NA)
      id = as.character(unseenList$animalNumber[i])
      change = append(CMRLOG[[id]]@stateList,unseenState)
      CMRLOG[[id]]@stateList <<- change
    }
  }
}

#State in the data is same as "not seen"
updateLogDead = function(t){
  for(i in 1:length(CMRLOG)){
    timesRecorded = length(CMRLOG[[i]]@stateList)
    timesNeedRecord = t - CMRLOG[[i]]@yearFirstCaptured + 1 
    timeDeficit = timesNeedRecord - timesRecorded
    if(timeDeficit == 1){
      unseenState = new("CMRState", time = t, wasSeen = FALSE, wasCaptured = FALSE,
                        currentSite = 0, isSeropositive = NA)
      change = append(CMRLOG[[i]]@stateList, unseenState)
      CMRLOG[[i]]@stateList <<- change
      }
    else if(timeDeficit != 0){
      a = CMRLOG[[i]]
      stop("something is wrong with recording deads and bird") 
    }
  }
}

#For time 0 CMR marking
markFirstTime = function(cmrParamMatrix, population){
 a = ncol(cmrParamMatrix)
 for(i in 1:ncol(cmrParamMatrix)){
   siteInds = getSitePop(i, population)
   markedInds = markNewIndividuals(siteInds, cmrParamMatrix[i,]$numberMarked, i, 0)
   h = population[is.element(population$animalNumber, markedInds$animalNumber),]
   population[is.element(population$animalNumber, markedInds$animalNumber),
              "isMarked"] = TRUE
   } 
 return(population)
}


#Takes a current population and returns the state that the population will be 
#at during the next time step
doTimeStep = function(totalPop, #type Population
                      epiParamMatrix, #data frame composed of type EpiSiteVec
                      currentTime, #numeric
                      dispersalMatrix){ #2D matrix of dispersal probabilities
  dispersedPop = c()
  for(i in 1:nrow(epiParamMatrix)){ #To clarify, there is one param row per site
    dispersedSite = dispersalTransitions(epiParamMatrix[i,], totalPop, i, 
                                  dispersalMatrix[i,])
    dispersedPop = rbind(dispersedPop, dispersedSite)
  }
  updatedPop = c()
  for(i in 1:nrow(epiParamMatrix)){
    updatedSite = postDispersalTransitions(epiParamMatrix[i,], dispersedPop, i)
    updatedPop = rbind(updatedPop, updatedSite)
  }
  return(updatedPop)
}


#Runs a simulation given parameters for sampling sites, runtime, and model input
runSim = function(popSizeVec, #atomic vector of nums. Length = number of sites.
                  #Entries = pop size at each site.
                  endTime, #numeric -- time steps to do
                  epiParamMatrix, #data frame composed of type EpiSiteVec
                  dispersalMatrix, #matrix holding dispersal probabilities
                  cmrParamMatrix){ #data frame composed of type MarkSiteVec
  population = newPopulation(popSizeVec)
  newestAnim <<- max(population$animalNumber)
  population = markFirstTime(cmrParamMatrix, population)
  simulationStats = getPopStats(population, length(popSizeVec), 0)
  for(i in 1:endTime){
    population = doTimeStep(population, epiParamMatrix, i, dispersalMatrix)
    population = captureMarkRecapture(population, cmrParamMatrix, i)
    timeStepStats = getPopStats(population, length(popSizeVec), i)
    simulationStats = rbind(simulationStats, timeStepStats)
  }
  #simCMRLog = CMRLOG
  #CMRLOG <<- list()
  
  return(simulationStats)
}


#UNIT TESTS... Admittedly, this code is not as polished as the rest.
#I didn't use a wide range of test cases, but code is simple. Should suffice.
#runTests = function(){
  #TEST POPULATION GENERATION
  sitePops = c(2, 5, 9, 30) #Four sites. Pop sizes of 2, 5, 9, 30.
  popGenerated = newPopulation(sitePops)
  #TEST GET STATS
  statsArr = getPopStats(popGenerated, 4, 1)
  compareVec = as.numeric(statsArr[,"populationSize"])
  stopifnot(identical(compareVec,sitePops))
  #TEST A TIME STEP
  epiVec1 = newEpiSiteVec(1, 2, .1, .1, .2)
  epiVec2 = newEpiSiteVec(2, 5, .1, .1, 0)
  epiVec3 = newEpiSiteVec(3, 9, .1, .1, .9)
  epiVec4 = newEpiSiteVec(4, 30, .1, .1, .5)
  disp = c(1/2, 1/2, 0, 0, #site 1
        1/5, 2/5, 1/5, 1/5, #site 2
        0, 1/9, 4/9, 4/9, #site 3
        0, 1/30, 4/30, 25/30) #site 4
  dispMatrix = matrix(disp, nrow = 4, ncol = 4, byrow = TRUE)
  epiMatrix = rbind(epiVec1, epiVec2, epiVec3, epiVec4)
  popStepped = doTimeStep(popGenerated, epiMatrix, 1, dispMatrix)
  #TEST DEATH
  epiVec = newEpiSiteVec(4, 30, .3, .1, .5) #reset to orig vals
  epiVec2 = newEpiSiteVec(1, 2, .3, .1, .2)
  sitePop = popGenerated[popGenerated$currentSite == 4,]
  sitePop = removeDying(epiVec, sitePop)
  #stopifnot(nrow(sitePop) == 21)
  #TEST DISPERSAL
  siteDispersalVec = c(.3,.3,.3,0)
  dispersed = disperse(sitePop, siteDispersalVec, epiVec)
  #stopifnot(nrow(dispersed[dispersed$currentSite == 2,]) == 6)
  #TEST CMR
  poppp = c(1,4,6)
  animNums = newPopulation(poppp)
  updateLogNew(animNums, 1, 1)
  seenNums = animNums
  updateLogSeen(seenNums,1, 1)
  getNum = as.character(1)
  print(CMRLOG[[getNum]]@stateList[[2]]@wasSeen == TRUE)
  #stopifnot(CMRLOG[[getNum]]@stateList[[2]]@wasSeen == TRUE)
  CMRLOG <<- list()
  #TESTSIM
  captureParam1 = newMarkSiteVec(1,1,.5,.5)
  captureParam2 = newMarkSiteVec(2,1,.5,.5)
  captureParam3 = newMarkSiteVec(3,1,.5,.5)
  captureParam4 = newMarkSiteVec(4,1,.5 ,.5)
  captureMat = rbind(captureParam1, captureParam2, captureParam3, captureParam4)
  output = runSim(sitePops, 40, epiMatrix, dispMatrix, captureMat)
#}
#runTests()
#DESIGN:
#create population
#run simulation
#do both cross-sectional and longitudinal sampling
#estimate and save params
#add functionality for sensitivity analysis but don't call it
