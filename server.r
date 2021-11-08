server <- function(input, output) {
  output$ICcelltext <- renderText({
    "The information content of sequence-period cells is displayed, where the information content of a cell
    is the variance of the treatment effect estimator if that cell is excluded divided by the variance if that 
    cell is included. There may be some cells which are entirely black:
    this implies that the model as specified cannot be fitted if this cell is excluded."
  })
  
  output$plot1<-renderPlot({
    
    Xdes <- SWdesmat(input$T)
    
    
    #If there are any periods with only 2 treatment sequences, with
    #differential exposure, cannot calculate the information content
    #of either of these two sequences.
    #Need to flag that the cells in these periods MUST be included and
    #thus do not have an information content.
    
    varmatall <- CRTVarGeneral(Xdes,input$m,input$rho0,input$r,input$type)
    varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
    
    for(i in 1:nrow(Xdes)){
      for (j in 1:ncol(Xdes)){
        if(is.na(Xdes[i,j])==TRUE)  varmat_excl[i,j] <- NA
        else if(is.na(Xdes[i,j])==FALSE) {
          Xdesij <- Xdes
          Xdesij[i,j] <- NA
          varmat_excl[i,j] <- CRTVarGeneral(Xdesij,input$m,input$rho0,input$r,input$type)/varmatall
        }
      }
    }
    
    varmat_excl<-round(varmat_excl, 4)
    melted_varmatexcl <- melt(varmat_excl)
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"
    
    color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))
    # if(sum(melted_varmatexcl$value==2.2772, na.rm=TRUE) > 0)    
    #   color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
    
    T <- ncol(Xdes)
    K <- nrow(Xdes)
    
    ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) +
      ggtitle("") +
      geom_tile(colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +
      geom_text(aes(Period, Sequence, label = round(value,2)), color = "black", size = 5) +
      scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)])
    
    
  })
}
