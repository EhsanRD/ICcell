---
  title: "Outputs"
output: html_document
runtime: shiny
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = '~/../')
```

## Embedded Shiny Application
```{r tabsets, echo=FALSE,message=FALSE}
shinyAppFile(appFile="appICcell/app.R",
             options = list(
               width = "100%", height = 700
             )
)
```

