library(shiny)
library(trackViewer)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)

colors <- c("blue", "orange", "green", "red", "purple", "brown", 
            "pink", "grey", "yellow", "cyan", "lightblue", "magenta", "violet")

#Data manipulation Part Starts
file <- file.choose()
mutationn1 <- read.delim(file, sep = "\t",header=TRUE)
mutationn1 <- as.data.frame(mutationn1)
mutationn_n <- mutationn1 %>% 
  separate_rows(Protein.change, sep = ",") %>%
  rename(amino_acid_change = Protein.change)

#Checking if NA values in column and clean it if occurs
dim1 <- nrow(mutationn_n)
if (dim1 != nrow(mutationn_n)) {
  warning("Some rows containing empty values in the 'amino_acid_change' column have been removed.")
}

mutationn_n <- mutationn_n[mutationn_n$amino_acid_change != "", ]
mutationn_n <- mutationn_n[, c("amino_acid_change","Clinical.significance..Last.reviewed.","Gene.s.")]
mutationn_n$Clinical.significance..Last.reviewed. <- gsub("(?=La).*", "", mutationn_n$Clinical.significance..Last.reviewed., perl = TRUE)
mutationn_n$Variant_Classification <- gsub("[(]", "", mutationn_n$Clinical.significance..Last.reviewed.)
mutationn_n$Variant_Classification <- gsub(" ", "_", mutationn_n$Variant_Classification)
mutationn_n$Gene.s.<- str_replace_all(mutationn_n$Gene.s., "[^[:alnum:]]", "")
mutationn_n$Hugo_Symbol <- gsub("LOC105371046", "", mutationn_n$Gene.s.)
mutationn_n <- mutationn_n[, c("amino_acid_change","Variant_Classification","Hugo_Symbol")]
mutationn_n$aaa <- as.numeric(gsub("[^[:digit:]]", "", mutationn_n$amino_acid_change))
#Data Manipulation part Ends

#User Interface Starts
ui <- fluidPage(
  titlePanel("Track Viewer"),
  
  sidebarLayout(
    
    #Sidebarpanel starts
    sidebarPanel(
      radioButtons("scs", "Select Clinical Significance",
                   choices = c("All", "Uncertain_significance", "Benign", "Conflicting_interpretations_of_pathogenicity","Likely_pathogenic","Likely_benign","Benign/Likely_benign","Pathogenic","risk_factor","Pathogenic/Likely_pathogenic")),
      
      actionButton("submit", "Submit"),
      
      sliderInput("range", "Select a range:",
                  min = 1, max = max(mutationn_n$aaa)+300,
                  value = c(1, max(mutationn_n$aaa)+300),
                  step = 1),
      
      h4("Choose Lollipop Colors"),
      selectInput("l1", "Uncertain_significance", choices = colors, selected = colors[1]),
      selectInput("l2", "Benign", choices = colors, selected = colors[2]),
      selectInput("l3", "Conflicting_interpretations_of_pathogenicity", choices = colors, selected = colors[3]),
      selectInput("l4", "Likely_pathogenic", choices = colors, selected = colors[4]),
      selectInput("l5", "Likely_benign", choices = colors, selected = colors[5]),
      selectInput("l6", "Benign/Likely_benign", choices = colors, selected = colors[6]),
      selectInput("l7", "Pathogenic", choices = colors, selected = colors[7]),
      selectInput("l8", "risk_factor:", choices = colors, selected = colors[8]),
      selectInput("l9", "Pathogenic/Likely_pathogenic", choices = colors, selected = colors[9])
      
    ),
    
    #Main Panel Starts
    mainPanel(
      h4(textOutput("count")),
      downloadButton("downloadSVG", "Download as SVG"),
      downloadButton("downloadPNG", "Download as PNG"),
      downloadButton("downloadJPEG", "Download as JPEG"),
      downloadButton("downloadTIFF", "Download as TIFF"),
      plotOutput("lollipopPlot")
    )
  )   
  
)
#User Interface Ends

#Server Starts
server <- function(input, output) {
  
  filter_data <- reactive({
    if (input$scs == "All") {
      return(mutationn_n)
    } else {
      return(mutationn_n[mutationn_n$Variant_Classification == input$scs,])
    }
  })
  
  output$count <- renderText({
    filtered_data <- filter_data()
    return(paste0(nrow(filtered_data), " amino acid change(s) found."))
  })
  
  #Submit Button Stars
  observeEvent(input$submit, {
    data_table <- filter_data()
    aaa <- data_table$aaa
    amino <- data_table$amino_acid_change
    gene_name <- unique(data_table$Hugo_Symbol)
    clinical <- unique(data_table$Variant_Classification)
    sample.gr <- GRanges("chr1", IRanges(aaa, width=1))
    features <- GRanges("chr1", IRanges(1, 1))
    legend <- list(labels=c("Uncertain_significance", "Benign", "Conflicting_interpretations_of_pathogenicity","Likely_pathogenic","Likely_benign","Benign/Likely_benign","Pathogenic","risk_factor","Pathogenic/Likely_pathogenic"), fill=c(input$l1,input$l2,input$l3,input$l4,input$l5,input$l6,input$l7,input$l8,input$l9))
    sample.gr$color <- ifelse(data_table$Variant_Classification == "Uncertain_significance", input$l1, ifelse(data_table$Variant_Classification == "Benign", input$l2, ifelse(data_table$Variant_Classification == "Conflicting_interpretations_of_pathogenicity",input$l3, ifelse(data_table$Variant_Classification=="Likely_pathogenic",input$l4, ifelse(data_table$Variant_Classification=="Likely_benign",input$l5, ifelse(data_table$Variant_Classification=="Benign/Likely_benign",input$l6,  ifelse(data_table$Variant_Classification=="Pathogenic",input$l7, ifelse(data_table$Variant_Classification=="risk_factor",input$l8, ifelse(data_table$Variant_Classification=="Pathogenic/Likely_pathogenic",input$l9, "white")))))))))
    sample.gr$SNPsideID <- sample(c("top","bottom"),size=length(aaa), replace=TRUE)
    names(sample.gr) <- paste0("p.",amino)
    
    
    #Outputs start
    output$lollipopPlot <- renderPlot({
      lolliplot(sample.gr,features,ylab=paste0(gene_name," ","Clinical Significance"),legend=legend,type="circle", ranges = GRanges("chr1", IRanges(input$range[1], input$range[2])))
      grid.text("ClinVar Color Scale of Clinical Significance Values",  
                x=0.5, y=0.9, 
                just="top",
                gp=gpar(cex=1.2, z=2))
    })
    
    output$downloadSVG <- downloadHandler(
      filename = function() {
        paste("lolliplot", Sys.Date(), ".svg", sep="")
      },
      content = function(file) {
        svg(filename = file,width = 1230, height = 434)
        lolliplot(sample.gr, features,ylab=paste0(gene_name," ","Clinical Significance"), legend=legend,type="circle", ranges = GRanges("chr1", IRanges(input$range[1], input$range[2])))
        grid.text("ClinVar Color Scale of Clinical Significance Values", x=.5, y=.9, just="top", 
                  gp=gpar(cex=1.2))
        dev.off()
      }
    )
    
    output$downloadPNG <- downloadHandler(
      filename = function() {
        paste("lolliplot", Sys.Date(), ".png", sep="")
      },
      content = function(file) {
        png(filename = file,width = 1230, height = 434)
        lolliplot(sample.gr, features,ylab=paste0(gene_name," ","Clinical Significance"), legend=legend,type="circle", ranges = GRanges("chr1", IRanges(input$range[1], input$range[2])))
        grid.text("ClinVar Color Scale of Clinical Significance Values", x=.5, y=.9, just="top", 
                  gp=gpar(cex=1.2))
        dev.off()
      }
    )
    
    output$downloadTIFF <- downloadHandler(
      filename = function() {
        paste("lolliplot", Sys.Date(), ".tiff", sep="")
      },
      content = function(file) {
        tiff(filename = file,width = 1230, height = 434)
        print(lolliplot(sample.gr, features,ylab=paste0(gene_name," ","Clinical Significance"), legend=legend,type="circle", ranges = GRanges("chr1", IRanges(input$range[1], input$range[2]))))
        grid.text("ClinVar Color Scale of Clinical Significance Values", x=.5, y=.9, just="top", 
                  gp=gpar(cex=1.2))
        dev.off()
      }
    )
    
    output$downloadJPEG <- downloadHandler(
      filename = function() {
        paste("lolliplot", Sys.Date(), ".jpeg", sep="")
      },
      content = function(file) {
        jpeg(filename = file,width = 1230, height = 434)
        lolliplot(sample.gr, features,ylab=paste0(gene_name," ","Clinical Significance"), legend=legend,type="circle", ranges = GRanges("chr1", IRanges(input$range[1], input$range[2])))
        grid.text("ClinVar Color Scale of Clinical Significance Values", x=.5, y=.9, just="top", 
                  gp=gpar(cex=1.2))
        dev.off()
      }
    )
    
  })
  #Submit button ends
}
#Server Ends

shinyApp(ui, server)
