library(shiny)
library(splines)

# https://statcomp2.app.vumc.org/rips/

source("predict.lm.R") # Maintain version it worked with (R 4.2.2)
load("cost.model.3.Rdata")
load("qaly.model.3.Rdata")
load("defaults.Rdata")



ages <- 22:60
defaults <- defaults[rep(1, length(ages)),]
defaults$IT <- ages

# Reference

noscreen <- defaults
noscreen$M <- rep(1, length(ages))

refcosts <- suppressWarnings(predict(cost.model.3, noscreen))
refqalys <- suppressWarnings(predict(qaly.model.3, noscreen))


# 
# vars <- read.csv("variables.csv")
# base <- as.numeric(vars$value)
# base <- data.frame(matrix(base, 1))
# excel <- c(outer(LETTERS, c('', LETTERS), function(a,b) paste0(b,a)))
# names(base) <- excel[2:(ncol(base)+1)] #vars$short
# 
# predict(cost.model.3, base)
#   
# 
# nas <- which(is.na(base))
# # Forward it's 241 where it starts working
# # Backward it's 103 where it starts working
# nas <- nas[nas >= 103 & nas <= 241]
# y <- base
# y[nas] <- base[1,nas]
# for(i in nas)
# {
#   y[i] <- NA
#   cat(i, "\n")
#   tryCatch(predict(cost.model.3, base), error=function(e) print(e))
#   y[i] <- base[1,i]
# }
# 
# 


ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("RISE Meta Model"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    tabsetPanel(
      tabPanel(
        "General",
        br(),
        sliderInput("sac", h3("Screening Assay Cost"), min=100, max=300, value=250), # Q
        sliderInput("cac", h3("Confirming Assay Cost"), min=100, max=300, value=186), # R
        sliderInput("ls.variant.prev", h3("Lynch Variant Prevalence (%)"), min=0.255, max=0.427, 0.341), # CO
        sliderInput("fh.variant.prev", h3("FH Variant Prevalence (%)"), min=0.235, max=0.611, 0.39), # EN
        sliderInput("bc.hboc.variant", h3("HBOC Variant Prevalence (%)"), min=0.404, max=0.603, 0.495), # V
        #        sliderInput("bc.prop.brac2", h3("Proportion of Carriers: BRCA2"), min=0.184, max=0.390, 0.289) # X
        
        # Utility Well was not run during the PSA
      ),
      tabPanel(
        "FH",
        br(),
        sliderInput("fh.fhp.190.220", h3("Relative Risk for MI FH(+), 190-220"), min=1.76, max=4.74, 2.86), # FE
        sliderInput("fh.notest.previouscve", h3("No Test, Previous CVE proportion MI"), min=0.47, max=0.783, 0.6), # IC
        sliderInput("fh.test.eff.size.20", h3("Test Effect Size, Age 20-40"), min=0.225, max=0.375, 0.3), #GN
      ),
      tabPanel(
        "HBOC",
        br(),
        # Breast Cancer 5yr survial was a multiple parameter calibration, most effect was column BC
        sliderInput("bc.mortal.rr.early", h3("BC Mortality Relative Reduction: Early Stage"), min=0.825, max=1.06, 0.943), # BC
        sliderInput("bc.utility.post", h3("Utility Post Breast Cancer"), min=0.72, max=0.9, 0.81),  # BI
        sliderInput("bc.hr.post.oo.brac1", h3("BC Hazard Ratio Post-Oophorectomy: BRCA1"), min=0.242, max=1.79, 0.63), # AT
        sliderInput("bc.hr.post.oo.brac2", h3("BC Hazard Ratio Post-Oophorectomy: BRCA2"), min=0.040, max=2.387, 0.36), # AU
        sliderInput("bc.int.screen.cost",  h3("Annual Intense Screening (MRI) Cost"), min=780, max=2065, 1403), # BQ
        sliderInput("bc.carriers.mri.prop", h3("Proportion of Known Carriers Who Get MRI"), min=0.37, max=0.96, 0.75),  # AH
      ),
      tabPanel(
        "Lynch",
        br(),
        wellPanel(
          sliderInput("ls.a.crc.surv", h3("Stage A CRC Surveilled"), min=0.078, max=0.930, 0.500),  # DE 
          sliderInput("ls.b.crc.surv", h3("Stage B CRC Surveilled"), min=0.043, max=0.86, 0.43), # DF
          sliderInput("ls.c.crc.surv", h3("Stage C CRC Surveilled"), min=0.043, max=0.86, 0.43), # DG
          sliderInput("ls.d.crc.surv", h3("Stage D CRC Surveilled"), min=0.043, max=0.86, 0.43) # DH
        ),
        wellPanel(
          sliderInput("ls.a.crc.nonsurv", h3("Stage A CRC Non-Surveilled"), min=0.09, max=0.66, 0.35), # DA
          sliderInput("ls.b.crc.nonsurv", h3("Stage B CRC Non-Surveilled"), min=0.09, max=0.66, 0.35), # DB
          sliderInput("ls.c.crc.nonsurv", h3("Stage C CRC Non-Surveilled"), min=0.09, max=0.66, 0.35), # DC
          sliderInput("ls.d.crc.nonsurv", h3("Stage D CRC Non-Surveilled"), min=0.017, max=0.507, 0.19) # DD
        ),
        sliderInput("ls.no.surv.crc.risk.70m", h3("No Surveillance CRC Risk by 70, Males"), min=0.087, max=0.865, 0.464), # CX
        sliderInput("ls.screen.cost", h3("Colonscopy Screening Cost/Yr"), min=755, max=2285, 1555), # DU
        sliderInput("ls.hr.surv.nosurv", h3("CRC Hazard Ratio Surveillance vs. None"), min=0.219, max=0.670, 0.387) # CY
      )
    )
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    h3("Results"),
    plotOutput("qalyPlot"),
    plotOutput("costPlot"),
    plotOutput("icerPlot"),
    p("This meta-model represents the model used in Guzauskas GF, et al., 'Population Genomic Screening for Three Common Hereditary Conditions : A Cost-Effectiveness Analysis'. Ann Intern Med. 2023 May;176(5):585-595. doi: 10.7326/M22-0846. It is a regression of the probability sensitivity analysis data created via block selection using splines and interactions in R. The code is published at",
      a(href="https://github.com/spgarbet/risesim", "https://github.com/spgarbet/risesim"), 
      "."),
    p("Built using ",
      a(href="https://www.r-project.org/", img(src="Rlogo.png", alt="R", width="60")),
      "with the ",
      a(href="https://cran.r-project.org/web/packages/shiny/shiny.pdf", "shiny"),
      "package."
    )
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output)
{
  v <- reactiveValues(data = defaults, costs=refcosts, qalys=refcosts, sd=10)
  
  observe({
    #v$data$P <- input$well
    v$data$Q  <- input$sac
    v$data$R  <- input$cac
     
    v$data$BI <- input$bc.utility.post
    v$data$BQ <- input$bc.int.screen.cost
    v$data$BC <- input$bc.mortal.rr.early
    v$data$V  <- input$bc.hboc.variant / 100
    v$data$AU <- input$bc.hr.post.oo.brac2
    v$data$AT <- input$bc.hr.post.oo.brac1
    v$data$AH <- input$bc.carriers.mri.prop

# FIXME
# Model explosion here
#    v$data$X  <- input$bc.prop.brac2
#    cat(input$bc.prop.brac2, "\n")
    
    v$data$EN <- input$fh.variant.prev /100
    v$data$FE <- input$fh.fhp.190.220
    v$data$IC <- input$fh.notest.previouscve
    v$data$GN <- input$fh.test.eff.size.20

    v$data$CO <- input$ls.variant.prev / 100
    v$data$CX <- input$ls.no.surv.crc.risk.70m
    v$data$DU <- input$ls.screen.cost
    v$data$CY <- input$ls.hr.surv.nosurv
    v$data$DE <- input$ls.a.crc.surv
    v$data$DD <- input$ls.d.crc.nonsurv
    v$data$DG <- input$ls.c.crc.surv
    v$data$DB <- input$ls.b.crc.nonsurv
    
    v$costs <- suppressWarnings(predict(cost.model.3, v$data))
    v$qalys <- suppressWarnings(predict(qaly.model.3, v$data))
  })

  output$qalyPlot <- renderPlot({
    plot(x=ages, y=100000*(v$qalys), type="p",
         xlab="Age at Time of Screening (y)", ylab="Incremental QALY per 100,000",
         ylim=c(100, 550))
  })
  
  output$costPlot <- renderPlot({
    plot(x=ages, y=v$costs/10, type="p",
         xlab="Age at Time of Screening (y)", ylab="Incremental Cost per 100,000 ($ millions)",
         ylim=c(15,45))
  })
  output$icerPlot <- renderPlot({
    #icer <- (v$costs-refcosts)/(v$qalys-refqalys)/10000
    icer <- (v$costs)/(v$qalys)/1000
#    ymin <- min(c(-1, icer[is.finite(icer)]))
#    ymax <- max(c(1,  icer[is.finite(icer)]))
    plot(x=ages, y=icer, type="p", xlab="Age (y)", ylab="Cost per QALY Gained ($ thousands)",
         ylim=c(0, 200))
    abline(h=100)
  })
}

shinyApp(ui, server)
