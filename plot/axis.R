drawsideaxis <- function(axisat = NULL, axislabels = NULL,
                         axisline = 0, axislwd = 1, axislas = 0,
                         axiscex = 1, labelpos = -0.5){
    

    
    if(!is.null(axisat)) {
        if(!is.null(axislabels)) {
            ## axis(4,
            ##      at = axisat,
            ##      labels = axislabels,
            ##      tick = FALSE,
            ##      xpd = NA,
            ##      line = labelpos,
            ##      las = axislas,
            ##      cex.axis = axiscex)
        } else {
            axis(4,
                 at = axisat,
                 ## labels = ylabels,
                 tick = FALSE,
                 xpd = NA,
                 line = labelpos,
                 las = axislas,
                 cex.axis = axiscex)
        }
        axis(2,
             at = axisat,
             labels = FALSE,
             tck = 0.05,
             line = axisline,
             lwd = axislwd)
        
        axis(4,
             at = axisat,
             labels = FALSE,
             tck = 0.05,
             line = axisline,
             lwd = axislwd)
    } else {
        axis(2,
             labels = FALSE,
             tck = 0.05,
             line = axisline,
             lwd = axislwd)
        
        axis(4,
             labels = FALSE,
             tck = 0.05,
             line = axisline,
             lwd = axislwd)

        ## axis(4,
        ##      tick = FALSE,
        ##      xpd = NA,
        ##      line = labelpos,
        ##      las = 2,
        ##      cex.axis = axiscex)        
    }
}
