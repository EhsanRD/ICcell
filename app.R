library(shiny)
library(ggplot2)
library(reshape2)
library(plyr)
library(swCRTdesign)
library(matrixcalc)
library(scales)
library(tidyverse)
library(shinythemes)
library(plotly)
library(Matrix)
library(reshape2)
library(ggplot2)

#setwd("~/Google Drive/Shared drives/Ehsan PhD work/Codes")
#source("ICcell_func.R", local=TRUE)

# Functions for generating design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}

#Designs missing cluster-period cells:

###this calculates the variance of the treatment effect estimator for the given inputs:###

#General function (taken from ExpDecayVar.r in Information content folder)
CRTVarGeneral <- function(Xmat, m, rho0, r, type) {

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP

  sig2 <- sig2E/m

  T <- ncol(Xmat)
  K <- nrow(Xmat)

  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])

  #Variance matrix for one cluster, with decay in correlation over time
  #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag(sig2 +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  #Variance matrix for all clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]

  #vartheta <- solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)]

  #xtry <- try(solve((t(Zmat)%*%solve(Vall)%*%Zmat))) 
  #if('try-error' %in% class(xtry)) return(NA)

  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])

  #return(vartheta)
}


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel(h1("",h2("Using Jessica Kasza GitHub codes & words"),h3("Information content of cluster-period cells of SW designs"))),
  # Sidebar with a slider input for number of bins
  
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "T", label = "Number of periods:",
                  min = 2, max = 20, value = 4,  step = 1),
      numericInput("m",
                   "Number of subjects in each cluster-period, m:",
                   min = 1,
                   max=1000,
                   step = 1,
                   value = 100),
      numericInput(
        "rho0", "Intra-cluster correlation",
        min = 0,
        max=0.2,
        step = 0.001,
        value = 0.05),
      radioButtons("type", label = ("Allow for decay correlation"), 
                   choices = list("Yes" = 1, "No" = 0), selected = 1),
      numericInput("r",
                   "Cluster auto-correlation",
                   min = 0,
                   max=1,
                   step = 0.05,
                   value = 0.95)),
  #   conditionalPanel(
  #     condition = "input.type == '1'",
  #     numericInput("r",
  #                  "Cluster auto-correlation",
  #                  min = 0,
  #                  max=1,
  #                  step = 0.05,
  #                  value = 1)),
  # ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Information content of cells",  
                 plotOutput("plot1"),
                 textOutput("ICcelltext")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
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


# Run the application
shinyApp(ui = ui, server = server)
