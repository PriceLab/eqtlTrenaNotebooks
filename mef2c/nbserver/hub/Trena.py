import zmq, json
import pandas as pd
from ipyTrenaViz import *
import time, os

class Trena:

    def __init__(self, genomeName):
       socketContext = zmq.Context();
       self.trenaServer = socketContext.socket(zmq.REQ)
       self.trenaServer.connect("tcp://trena:%s" % "5548")
       self.tv = ipyTrenaViz()
       display(self.tv)
       self.tv.setGenome(genomeName)

    def version(self):
       return(1.02)

    def display(self):
       display(self.tv)

    def ping(self):
        print("sending ping to treanServer")
        msg = {'cmd': 'ping', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        return(response)

    def showGenomicRegion(self, regionString):
        self.tv.showGenomicRegion(regionString);

    def getGenomicRegion(self):
        return(self.tv.getBrowserState()["chromLocString"]);

    def dataFrameFrom3partList(self, list):
        data = list['tbl']
        rownames = list['rownames']
        colnames = list['colnames']
        df = pd.DataFrame(data)
        df.columns = colnames
        rownameList = {}
        for i in range(len(rownames)):
          rownameList[i] = rownames[i]
        df = df.rename(rownameList)
        return(df)

    def getBounds(self, asString=False):
        msg = {'cmd': 'getGenomicBounds', 'status': 'request', 'callback': '', 'payload': asString}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload)

    def getModelNames(self):
        msg = {'cmd': 'getModelNames', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload)

    def getModel(self, name):
        msg = {'cmd': 'getModel', 'status': 'request', 'callback': '', 'payload': name}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        if(response["status"] != "success"):
           print("no model found with name '%s'" % name)
           return(payload)
        payload = response["payload"]
        tblAsList = payload["tbl"]
        tbl = self.dataFrameFrom3partList(tblAsList)
        return(tbl)

    def getVariants(self, minScore, display, color, trackHeight=50):
        payload = {"roi": self.getGenomicRegion(), 'minScore': minScore}
        msg = {'cmd': 'getVariants', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        tblAsList = payload["tbl"]
        tbl = self.dataFrameFrom3partList(tblAsList)
        tbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(tbl, "variants >= %4.2f" % minScore, "SQUISHED", color, trackHeight)
        return(tbl)

    def findVariantsInModel(self, motifTrack, variantsTrack, candidateTFs, tfMotifMappingName, shoulder):
        roi = self.getGenomicRegion()
        payload = {"roi": roi,
                   'motifTrack': motifTrack,
                   'variantsTrack': variantsTrack,
                   'candidateTFs': candidateTFs,
                   'tfMotifMappingName': tfMotifMappingName,
                   'shoulder': shoulder}
        msg = {'cmd': 'findVariantsInModelForRegion', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        if(response["status"] != "success"):
           print("no variants in model in region %s" % roi)
           return(payload)
        payload = response["payload"]
        tblAsList = payload["tbl"]
        tbl = self.dataFrameFrom3partList(tblAsList)
        return(tbl)

    def getWholeGenomeVariantsInRegion(self, altToRefRatio, minAltCount, display, color, trackHeight=50):
        roi = self.getGenomicRegion()
        payload = {"roi": roi,
                   'altToRefRatio': altToRefRatio,
                   'minAltCount': minAltCount}
        msg = {'cmd': 'getWholeGenomeVariants', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        if(response["status"] != "success"):
           print("no variants found in region %s" % roi)
           return(payload)
        payload = response["payload"]
        tblAsList = payload["tbl"]
        tbl = self.dataFrameFrom3partList(tblAsList)
        tbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(tbl, "wgs variants >= %4.1f, %4d" % (altToRefRatio, minAltCount),
                                            "SQUISHED", color, trackHeight)
        return(tbl)

    def getExpressionMatrixNames(self):
        msg = {'cmd': 'getExpressionMatrixNames', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload)

    def summarizeExpressionMatrices(self):
        msg = {'cmd': 'summarizeExpressionMatrices', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(self.dataFrameFrom3partList(payload))

    def getFootprintsInRegion(self, display):
        payload = {"roi": self.getGenomicRegion()}
        msg = {'cmd': 'getFootprintsInRegion', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]

        tblAsList = payload["tbl"]
        regTbl = self.dataFrameFrom3partList(tblAsList)
        regTbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(regTbl, "footprints", "SQUISHED", "blue")
        return(regTbl)

    def getDHSinRegion(self, display):
        payload = {"roi": self.getGenomicRegion()}
        msg = {'cmd': 'getDHSRegions', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        if(response["status"] != "success"):
           print("no DHS regions found")
           return(payload)
        payload = response["payload"]
        tblAsList = payload["tbl"]
        regTbl = self.dataFrameFrom3partList(tblAsList)
        regTbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(regTbl, "DHS", "SQUISHED", "darkreen", trackHeight=50)
        return(regTbl)

    def getEnhancersInRegion(self, display):
        payload = {"roi": self.getGenomicRegion()}
        msg = {'cmd': 'getEnhancers', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        if(response["status"] != "success"):
           print("no enhancers found")
           return(payload)
        tblAsList = payload["tbl"]
        regTbl = self.dataFrameFrom3partList(tblAsList)
        regTbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(regTbl, "Enhancers", "SQUISHED", "darkgreen", trackHeight=50)
        return(regTbl)

    def getMotifsInRegion(self, motif, matchScore, display, color="blue"):
        payload = {"roi": self.getGenomicRegion(),
                   "motifs": motif,
                   "matchScore": matchScore}
        msg = {'cmd': 'getMotifsInRegion', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        if(response["status"] != "success"):
           print("no motifs found in %s at match >= %d" % (self.getGenomicRegion(), matchScore))
           return()
        tblAsList = payload["tbl"]
        regTbl = self.dataFrameFrom3partList(tblAsList)
        if(display):
           tbl = regTbl[[0,1,2,3,4]]
           self.tv.addBedTrackFromDataFrame(tbl, motif, "SQUISHED", color, trackHeight=50)

        return(regTbl)

    def getDHSMotifsinRegion(self, display):
        payload = {"roi": self.getGenomicRegion()}
        msg = {'cmd': 'getDHSMotifsInRegion', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]

        tblAsList = payload["tbl"]
        regTbl = self.dataFrameFrom3partList(tblAsList)
        regTbl.key = payload["key"]
        if(display):
           self.tv.addBedTrackFromDataFrame(regTbl, "DHS motifs", "SQUISHED", "magenta", trackHeight=50)
        return(regTbl)

    #def findVariantsInModel(self, modelName, shoulder, display):
    #    payload = {"modelName": modelName, "shoulder": shoulder};
    #    msg = {'cmd': 'findVariantsInModel', 'status': 'request', 'callback': '', 'payload': payload}
    #    self.trenaServer.send_string(json.dumps(msg))
    #    response = json.loads(self.trenaServer.recv_string())
    #    payload = response["payload"]
    #    tblAsList = payload["tbl"]
    #    varTbl = self.dataFrameFrom3partList(tblAsList)
    #    varTbl.key = payload["key"]
    #    if(display):
    #       self.tv.addBedTrackFromDataFrame(varTbl.loc[:, ['chrom', 'pos', 'pos', 'rsid']], "voi", "SQUISHED", "darkred")
    #    return(varTbl)

    def displayFootprints(self, url):
        self.tv.addBedTrackFromDataFrame(url)

    def addBedTrackFromDataFrame(self, tbl, trackName, trackMode, color, trackHeight=200):
        return(self.tv.addBedTrackFromDataFrame(tbl, trackName, trackMode, color, trackHeight))

    def displayGraphFromFile(self, filename, modelNames):
        self.tv.displayGraphFromFile(filename, modelNames)

    def setStyle(self, filename):
        self.tv.setStyle(filename)

    def sessionInfo(self):
        msg = {'cmd': 'getSessionInfo', 'status': 'request', 'callback': '', 'payload': ""}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload);

    def listSharedData(self):
        msg = {'cmd': 'listSharedData', 'status': 'request', 'callback': '', 'payload': ""}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload);

    def createGeneModel(self, targetGene,  solverNames, tbl_regRegions, tfMap, matrixName):

        payload = {'targetGene': targetGene,
                   'solverNames': solverNames,
                   'tblRegulatoryRegionsCacheKey': tbl_regRegions.key,   # used to look up in cache
                   'tfMap': tfMap,
                   'matrixName': matrixName}
        msg = {'cmd': 'createGeneModel', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        tblAsList = payload["tbl"]
        tbl = self.dataFrameFrom3partList(tblAsList)
        tbl.key = payload["key"]
        return(tbl)

    def displayMultiModelGraph(self, targetGene, modelList):
       modelNames = list(modelList.keys())
       for modelName in modelNames:
          print(' now reducing modelName %s' % modelName)
          tbl = modelList[modelName]['model']
          modelList[modelName]['model'] = tbl.key
          tbl = modelList[modelName]['regions']
          modelList[modelName]['regions'] = tbl.key

       payload = {"targetGene": targetGene, "models": modelList};
       msg = {'cmd': 'buildMultiModelGraph', 'status': 'request', 'callback': '', 'payload': payload}
       self.trenaServer.send_string(json.dumps(msg))
       response = json.loads(self.trenaServer.recv_string())
       g_json = response["payload"]
       open("g.json", "w")
       f = open("g.json", "w")
       f.write(g_json)
       f.close()
       self.displayGraphFromFile("g.json", modelNames)
       print("after calling displayGraphFromFile");
       #return(payload)

    def createTaggedDataFrame(self, rows, columns):
        payload = {'rows': rows, 'cols': columns}
        msg = {'cmd': 'createTaggedDataFrame', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        tblAsList = payload["tbl"]
        pTbl = self.dataFrameFrom3partList(tblAsList)
        pTbl.key = payload["key"]
        return(pTbl)

    def findTaggedDataFrameOnServer(self, tbl):
        payload = tbl.key
        msg = {'cmd': 'identifyTaggedDataFrame', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload)

    def setWidgetHeight(self, newHeight):
        self.tv.setWidgetHeight(newHeight)


