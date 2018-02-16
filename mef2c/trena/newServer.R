library(MEF2C.data)
library(SingleGeneAnalyzer)
library(rzmq)
library(jsonlite)
library(trena)
library(trenaViz) # used here only for buildMultiModelGraph and addGeneLayout
library(RUnit)
PORT <- 5548
sessionInfo()
#------------------------------------------------------------------------------------------------------------------------
mef2c.data <- MEF2C.data()
# todo: make tss an argument to MEF2C.data ctor.  targetGene is provided by MEF2C.data::getTargetGene
sga <- SingleGeneAnalyzer(genomeName="hg38", targetGene="MEF2C", targetGene.TSS=88904257, mef2c.data)
cache <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
handleMessage <- function(msg)
{
   printf("incoming message '%s' at %s", msg$cmd, date())

   if(msg$cmd == "summarizeExpressionMatrices"){
      tbl.summary <- summarizeExpressionMatrices(sga)
      tbl.summary.as.list <- dataFrameToPandasFriendlyList(tbl.summary)
      response <- list(cmd=msg$callback, status="success", callback="", payload=tbl.summary.as.list)
      }
   else if(msg$cmd == "getGenomicBounds"){
      print("dispatching getGenomicBounds, msg: ")
      #print(msg)
      asString = msg$payload
      result <- getGenomicBounds(mef2c.data, asString)
      response <- list(cmd=msg$callback, status="success", callback="", payload=result);
      }
   else if(msg$cmd == "getFootprintsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roiString <- msg$payload$roi
      tbl.reg <- SingleGeneAnalyzer::getFootprintsForRegion(sga, roiString)
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.reg
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "getMotifsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      stopifnot("motifs" %in% names(msg$payload))
      stopifnot("matchScore" %in% names(msg$payload))
      roiString <- msg$payload$roi
      motifs <- msg$payload$motifs
      matchScore <- msg$payload$matchScore
      tbl.reg <- SingleGeneAnalyzer::findMotifsInRegion(sga, roiString, motifs, matchScore)
      printf("found %d motifs in sequence: %s", nrow(tbl.reg), roiString)
      if(nrow(tbl.reg) == 0){
         response <- list(cmd=msg$callback, status="failure", callback="", payload="no motifs found");
         }
      else{
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.reg
         tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
         payload <- list(tbl=tbl.fp.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         }
      } # getMotifsInRegion

   else if(msg$cmd == "getVariants"){
      roi.string <- msg$payload$roi
      source.name <- msg$payload$source.name
      tracking.name <- msg$payload$tracking.name
      supported.source.names <- c("IGAP.snpChip", "ADNI.WGS", "MAYO.eqtl.snps")
      if(!source.name %in% c(supported.source.names)){
         error.msg <- sprintf("source.name '%s' for variants unrecognized, should be one of %s",
                              paste(supported.source.names, collapse=", "))
         print(error.msg)
         response <- list(cmd=msg$callback, status="failure", callback="", payload=error.msg)
         }
      else{
         score.1.threshold <- msg$payload$score.1.threshold
         score.2.threshold <- msg$payload$score.2.threshold
         score.3.threshold <- msg$payload$score.3.threshold
         tbl.var <- getVariantsForRegion(sga, source.name, tracking.name, roi.string,
                                         score.1.threshold, score.2.threshold, score.3.threshold)
         if(nrow(tbl.var) == 0){
            response <- list(cmd=msg$callback, status="failure", callback="", payload="no variants in region");
            }
         else{
           printf("--- colnames of tbl.var: %s", paste(colnames(tbl.var), collapse=","))
           rownames(tbl.var) <- NULL
           key <- as.character(as.numeric(Sys.time()) * 100000)
           tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
           payload <- list(tbl=tbl.fp.as.list, key=key)
           response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
            } # nrow(tbl.snp) > 0
         } # else: good source.name provided
      } # getVariants

   else if(msg$cmd == "getMotifs"){
      roi.string <- msg$payload$roi
      source.name <- msg$payload$source.name
      tracking.name <- msg$payload$tracking.name
      score.threshold <- msg$payload$score.threshold
      supported.source.names <- c("allDNA-jaspar2018-human-mouse-motifs")
      if(!source.name %in% c(supported.source.names)){
         error.msg <- sprintf("source.name '%s' for motifs unrecognized, should be one of %s",
                              paste(supported.source.names, collapse=", "))
         print(error.msg)
         response <- list(cmd=msg$callback, status="failure", callback="", payload=error.msg)
         }
      else{
         tbl.motifs <- getMotifsForRegion(sga, source.name, tracking.name, roi.string, score.threshold)
         if(nrow(tbl.motifs) == 0){
            response <- list(cmd=msg$callback, status="failure", callback="", payload="no motifs in region");
            }
         else{
           printf("--- colnames of tbl.motifs: %s", paste(colnames(tbl.motifs), collapse=","))
           rownames(tbl.motifs) <- NULL
           key <- as.character(as.numeric(Sys.time()) * 100000)
           tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.motifs)
           payload <- list(tbl=tbl.fp.as.list, key=key)
           response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
            } # nrow(tbl.snp) > 0
         } # else: good source.name provided
      } # getMotifs

   else if(msg$cmd == "intersectTracks"){
      roiString   <- msg$payload$roi
      trackName.1 <- msg$payload$trackName_1
      trackName.2 <- msg$payload$trackName_2
      shoulder    <- msg$payload$shoulder
      tbl.intersect <- intersectTracks(sga, roiString, trackName.1, trackName.2, shoulder)
      if(nrow(tbl.intersect) == 0){
         errorMessage <- sprintf("no intersection of %s and %s in region %s", trackName.1, trackName.2, roiString);
         response <- list(cmd=msg$callback, status="failure", callback="", payload=errorMessage)
         }
      else{
        key <- as.character(as.numeric(Sys.time()) * 100000)
        tbl.as.list <- dataFrameToPandasFriendlyList(tbl.intersect)
        payload <- list(tbl=tbl.as.list, key=key)
        response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
        } # nrow(tbl.variants) > 0
      } # intersectTracks

   else if(msg$cmd == "getWholeGenomeVariants"){
      roiString <- msg$payload$roi
      altToRefRatio  <- msg$payload$altToRefRatio
      minAltCount <- msg$payload$minAltCount
      tbl.variants <- getWholeGenomeVariantsForRegion(sga, roiString, altToRefRatio, minAltCount)
      printf("whole genome tbl.variants: %d x %d", nrow(tbl.variants), ncol(tbl.variants))
      if(nrow(tbl.variants) == 0){
         response <- list(cmd=msg$callback, status="failure", callback="", payload="no variants in region");
         }
      else{
        key <- as.character(as.numeric(Sys.time()) * 100000)
        tbl.variants.as.list <- dataFrameToPandasFriendlyList(tbl.variants)
        payload <- list(tbl=tbl.variants.as.list, key=key)
        response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
        } # nrow(tbl.variants) > 0
      } # getVariants
   else if(msg$cmd == "getDHSRegions"){
      stopifnot("roi" %in% names(msg$payload))
      thresholdScore <- msg$payload$minScore
      roi = msg$payload$roi
      tbl.dhs.roi <- getDHSForRegion(sga, roi)
      printf("dhs regions found: %d", nrow(tbl.dhs.roi))
      response <- list(cmd=msg$callback, status="failure", callback="", payload="no overlapping regions");
      if(nrow(tbl.dhs.roi) > 0){
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.dhs.roi
         tbl.dhs.as.list <- dataFrameToPandasFriendlyList(tbl.dhs.roi)
         payload <- list(tbl=tbl.dhs.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         } # else: some overlap found
      } # getDHsRegions
   else if(msg$cmd == "getEnhancers"){
      stopifnot("roi" %in% names(msg$payload))
      thresholdScore <- msg$payload$minScore
      roi = msg$payload$roi
      tbl.dhs.roi <- getEnhancersForRegion(sga, roi)
      printf("dhs regions found: %d", nrow(tbl.dhs.roi))
      response <- list(cmd=msg$callback, status="failure", callback="", payload="no overlapping regions");
      if(nrow(tbl.dhs.roi) > 0){
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.dhs.roi
         tbl.dhs.as.list <- dataFrameToPandasFriendlyList(tbl.dhs.roi)
         payload <- list(tbl=tbl.dhs.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         } # else: some overlap found
      } # getDHsRegions
   else if(msg$cmd == "getSessionInfo"){
      info <- as.character(sessionInfo())
      response <- list(cmd=msg$callback, status="success", callback="", payload=info)
      }
   else if(msg$cmd == "getExpressionMatrixNames"){
      info <- sort(names(expression.matrix.files))
      response <- list(cmd=msg$callback, status="success", callback="", payload=info)
      }
   else if(msg$cmd == "getModelNames"){
      model.names <- getRegulatoryModelNames(sga)
      printf("model.names: %s", paste(model.names, collapse=", "))
      response <- list(cmd=msg$callback, status="success", callback="", payload=model.names)
      }
   else if(msg$cmd == "getModel"){
      modelName <- msg$payload
      tbl <- getRegulatoryModel(sga, modelName)
      response <- list(cmd=msg$callback, status="failure", callback="",
                       payload=sprintf("model '%s' not found", modelName))
      if(nrow(tbl) > 0){
         tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl)
         payload <- list(tbl=tbl.fp.as.list)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         } # if rows in tbl
      } # getModel
   else if(msg$cmd == "getDHSMotifsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roi <- msg$payload$roi
      tbl.reg <- getDHSMotifs(roi)
      if(nrow(tbl.reg) == 0){
         response <- list(cmd=msg$callback, status="failure", callback="", payload="no motifs in region");
         }
      else{
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.reg
         tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
         payload <- list(tbl=tbl.fp.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         }
      } # getDHSMotifsInRegion
   else if(msg$cmd == "findVariantsInModelForRegion"){
       printf("executing **findVariantsInModelForRegion**")
       stopifnot("roi" %in% names(msg$payload))
       stopifnot("motifTrack" %in% names(msg$payload))
       stopifnot("variantsTrack" %in% names(msg$payload))
       stopifnot("candidateTFs" %in% names(msg$payload))
       stopifnot("tfMotifMappingName" %in% names(msg$payload))
       stopifnot("shoulder" %in% names(msg$payload))
       response <- list(cmd=msg$callback, status="failure", callback="", payload="no variants in tf motifs in region");
       printf("about to call sga method")
       with(msg$payload, {
         printf("--- roi: %s", roi)
         printf("    motif.track: %s", motifTrack)
         printf("    variants.track: %s", variantsTrack)
         printf("    tfs: %s", paste(candidateTFs, collapse=", "))
         printf("    mapper: %s", tfMotifMappingName)
         printf("    shoulder: %d", shoulder)
         })
       tbl.var <- with(msg$payload, findVariantsInModelForRegion(sga,
                                              roi.string=roi,
                                              motif.track=motifTrack,
                                              variants.track=variantsTrack,
                                              candidate.tfs=candidateTFs,
                                              tfMotifMappingName=tfMotifMappingName,
                                              shoulder=shoulder))
      printf("findVariantsInModelForRegion, row count: %d", nrow(tbl.var))
      if(nrow(tbl.var) > 0){
         tbl.as.list <- dataFrameToPandasFriendlyList(tbl.var)
         payload <- list(tbl=tbl.as.list)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         } # if success
      } # findVariantsInModelForRegion
   else if(msg$cmd == "findVariantsInModel"){
      stopifnot("modelName" %in% names(msg$payload))
      stopifnot("shoulder" %in% names(msg$payload))
      modelName <- msg$payload$modelName
      shoulder <- msg$payload$shoulder
      tbl.var <- findVariantsInModel(modelName, shoulder, "footprints")
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.var
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "listSharedData"){
      filenames <- list.files("/home/trena/sharedData")
      response <- list(cmd=msg$callback, status="success", callback="", payload=filenames);
      }
   else if(msg$cmd == "createTaggedDataFrame"){
      rows <- msg$payload$rows
      cols <- msg$payload$cols
      tbl <- as.data.frame(matrix(runif(32), nrow=rows, ncol=cols, dimnames=list(LETTERS[1:rows],  letters[1:cols])))
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl
      tbl.for.pandas <- dataFrameToPandasFriendlyList(tbl)
      payload <- list(tbl=tbl.for.pandas, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "identifyTaggedDataFrame"){
      key <- msg$payload
      found <- key %in% names(cache)
      tbl.retrieved <- cache[[key]]
      printf("retrieved tbl with dimensions %d, %d", nrow(tbl.retrieved), ncol(tbl.retrieved))
      #print(tbl.retrieved)
      sum <- sum(as.matrix(tbl.retrieved))
      payload <- list(found=found, sum=sum)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "createGeneModel"){
      trena <- Trena("hg38")   # probably should create a single global instead
      payload <- msg$payload
      solver.names <- payload$solverNames
      printf("---- createGeneModel about to extract tbl.motifs from cache")
      key <- payload$tblRegulatoryRegionsCacheKey
      #printf("    key: %s", key)
      found <- key %in% names(cache)
      #printf("    found? %s", found)
      tbl.motifs <- cache[[key]]
      #printf("    tbl.motifs: %d, %d", nrow(tbl.motifs), ncol(tbl.motifs))
      tfMap <- payload$tfMap
      stopifnot(all(c("targetGene", "matrixName") %in% names(msg$payload)))
      targetGene <- toupper(msg$payload$targetGene)
      matrixName <- msg$payload$matrixName
      #tbl.motifs <- read.table("/home/trena/sharedData/tbl.bed", sep="\t", as.is=TRUE, stringsAsFactors=FALSE)
      tbl.motifs <- tbl.motifs[, 1:5]
      colnames(tbl.motifs) <- c("chrom", "start", "end", "motifName")
      #print(head(tbl.motifs))
      tbl.motifs$motifName <- sub("_well_16", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_well_20", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_hint_16", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_hint_20", "", tbl.motifs$motifName)
      #print(head(tbl.motifs))

      tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source=tfMap, expand.rows=FALSE)
      #print(head(tbl.motifs.tfs))

      mtx <- expression.matrices[[matrixName]]
      stopifnot(targetGene %in% rownames(mtx))
      print(dim(mtx))
      print(mean(mtx[targetGene,]))
      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs.tfs, mtx)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.geneModel
      print(tbl.geneModel)
      payload=list(tbl=dataFrameToPandasFriendlyList(tbl.geneModel), key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload)
      } # createGeneModel
   else if(msg$cmd == "buildMultiModelGraph"){
      printf("----- buildMultiModelGraph")
      payload = msg$payload
      print(payload)
      targetGene <- payload$targetGene
      models <- payload$models
      models.decached <- list()
      for(name in names(models)){
         model <- models[[name]]
         tbl.geneModel <-  cache[[model$model]]
         tbl.regions <- cache[[model$regions]]
         tbl.regions <- subset(tbl.regions, geneSymbol %in% tbl.geneModel$gene)
         models.decached[[name]]$model <- tbl.geneModel
         models.decached[[name]]$regions <- tbl.regions
         } # for model
      printf("---- after decaching")
      print(models.decached)
      g <- trenaViz::buildMultiModelGraph(targetGene, models.decached)
      xCoordinate.span <- 1500
      g.lo <- trenaViz::addGeneModelLayout(g, xPos.span=xCoordinate.span)
      printf("---- g.lo")
      print(g.lo)
      g.json <- trenaViz:::.graphToJSON(g.lo)
      response <- list(cmd=msg$callback, status="success", callback="", payload=g.json)
      } # buildMultiModelGraph
   else{
      response <- list(cmd=msg$callback, status="success", callback="",
                       payload=sprintf("unrecognized command: '%s'", msg$cmd))
      }

   print("--- about to leave handleMessage, response:")
   #print(response)

   response

} # handleMessage
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) {

   context = init.context()
   socket = init.socket(context,"ZMQ_REP")
   bind.socket(socket, sprintf("tcp://*:%d", PORT))

   errorFunction <- function(condition){
     printf("==== exception caught ===")
     print(as.character(condition))
     response <- list(cmd="handleError", status="error", callback="", payload=as.character(condition));
     send.raw.string(socket, toJSON(response))
     };

   while(TRUE) {
      tryCatch({
        printf("top of receive/send loop")
        raw.message <- receive.string(socket)
        msg = fromJSON(raw.message)
        printf("--- msg:")
        print(msg)

        stopifnot(is.list(msg))
        expected.fields <- c("cmd", "status", "callback", "payload")
        stopifnot(all(expected.fields %in% names(msg)))

        response <- handleMessage(msg)
        json.string <- toJSON(response, dataframe="values")
        send.raw.string(socket, json.string)
        Sys.sleep(1)
        }, error=errorFunction); # tryCatch
     } # while (TRUE)

} # if !interactive()
#------------------------------------------------------------------------------------------------------------------------
