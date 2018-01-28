.MEF2C.data <- setClass ("MEF2C.data",
                            representation = representation(
                               quiet="logical"
                               )
                            )

PORT <- 5548
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getExpressionMatrices', signature='obj', function(obj) standardGeneric ('getExpressionMatrices'))
#------------------------------------------------------------------------------------------------------------------------
MEF2C.data = function(quiet=TRUE)
{
   obj <- .MEF2C.data(quiet=quiet)
   obj

} # constructor
#----------------------------------------------------------------------------------------------------
setMethod('getExpressionMatrices', 'MEF2C.data',

    function(obj){
       mtx.names <- load(system.file(package="MEF2C.data", "extData", "mtx.withDimers.cer.ros.tcx.RData"))
       expression.matrices <- list()
       for(matrix.name in mtx.names){
          mtx <- eval(parse(text=matrix.name))
          expression.matrices[[matrix.name]] <- mtx
          }
       invisible(expression.matrices)
       })

#----------------------------------------------------------------------------------------------------
