# Exploring Transcriptional Regulation of the FRD3 gene in Athaliana, using trena
# We start by:
#
#   - loading a few R packages
#   - initializing some important variables
#   - launching trenaViz - which provides an igv genome browser
#     and cytoscape network visualization.
#
# Note that trenaViz appears in another web browser window, separate
# from this RStudio notebook window.  For most people using most
# browsers, trenaViz appears in a new browser tab.  You may wish to
# "tear off" and place it side-by-side with the RStudio window.  This
# allows you to see and execute the simple code which follows below
# AND to see and manipulate the interactive graphics provided by
# trenaVz.
#
# igv is useful for setting up parameters for your model building:
# i.e., setting the chromosomal region where transcription factor
# binding sites will be found.  cytoscape displays the resulting
# regulatory models - possibly more than one at at time for easy
# comparison of the consequences of your different modeling
# assumptions, as we will see below.
#
# Not_: If your web browser is configured to block popup windows, then
# the new trenaViz window/browser tab will not appear.  Each browser
# provides some means to allow popups - please enable that if you need
# to.
#

library(trena)
library(trenaViz)
library(FRD3.data)

PORTS <- 10000:10020   # trenaViz claims a websocket port in this range, over which R, igv and cytoscape communicate

tv <- trenaViz(portRange=PORTS)
setGenome(tv, "tair10")
frd3 <- FRD3.data()
target.orf <- frd3@misc.data$targetGene
setBrowserWindowTitle(tv, "FRD3")

# Among the many items provided in the _FR3D.data_ package loaded above
# is our current guess of the chromosomal bounds of a region of interest
# ("roi") for the study of this gene.  Get that value, then convert to
# the sort of "chrom:start-end" format igv understands, and ask igv (via
# trenaViz, *tv*) to display it.


roi <- getGenomicBounds(frd3)
roi.string <- with(roi, sprintf("%s:%d-%d", chrom, start, end))
showGenomicRegion(tv, roi.string)   # about 3kb upstream and downstream of primary tss


# Look for, then click the "minus" icon at the top right of the igv
# window.  This will zoom out by a factor of two, providing a regional
# sense of the distribution of genes and DHS region in and around FRD3.
# This view is worth keeping in mind when you pick a region for
# modeling.

# Return to the initial roi, 5875 bp more or less centered on the FRD3 primary TSS.
#
# Create a first model using a very permissive DHS threshold: any
# region with a score greater than 0.5 is evaluated to see if contains
# an athaliana transcription factor binding site.

m1 <- makeModelForRegion(frd3, dhs.cutoff=1.0, region=roi.string, trenaViz=tv)

# Examine the two data.frames returned: one (the model) describes the
# TF/gene relationships.  The other ("regions") describes the putative
# transcription factor binding sites in the DHS regions.  TFs mapped to
# those possible binding sites are *candidate regulators*.  trena
# evaluates their gene expression, and that of the target gene, to
# predict regulatory relationshiops.  We can see that *ARR10* is a
# strong contender.

View(m1$regions)

View(m1$model)

# Now make another model, this one (at least speculatively) specific to
# the less common transcript, visible as the second transcript in the
# FRD3 gene model displayed in igv.


roi.2 <- "3:2,571,751-2,573,886"
showGenomicRegion(tv, roi.2)
m2 <- makeModelForRegion(frd3, dhs.cutoff=1.0, region=roi.2, trenaViz=tv)

# Notice that the two models are similar, though the TFs are in a somewhat different order.

View(m2$model)

# Let's now look at these two models in cytoscape.  The trenaViz view will shift from "IGV" to "TRN".  Use the "Cycle Models" button, displayed above the trn network, to flip between the two models.

models <- list(left=m1, right=m2)
g <- buildMultiModelGraph(target.orf, models)
g.lo <- addGeneModelLayout(g, xPos.span=1500)
setGraph(tv, g.lo, names(models))
setStyle(tv, system.file(package="FRD3.data", "extdata", "style.js"))

# add a genome track (in the IGV tab) for all the occurences of 
# 
tbl.model <- m1$model
tbl.motifs <- m1$regions
# choose a tf with multiple binding sites
max.bindingSites <- max(tbl.model$bindingSites)
tf.with.max.bindingSites <- which(tbl.model$bindingSites == max.bindingSites)[1]
tf <- tbl.model$gene[tf.with.max.bindingSites]
tbl.bed <- motifTrackForTF(frd3, tbl.motifs, tf, tv)
