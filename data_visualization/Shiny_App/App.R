# This Shiny app is build to illustrate sample Toxic R visualization
# Author: Sooyeong Lim
# Date: 2020/12/07

# Required packages
library(tidyverse)
library(ggplot2)
library(shiny)
library(plotly)
library(scales)
library(ToxicR)

# Load functions for dichotomous and continous cases
# source("data_visualization/Shiny_App/dichotomous_functions.R")
# source("data_visualization/Shiny_App/continous_functions.R")
# 
# # Load plot generators
# source("data_visualization/Shiny_App/plot_dichotomous_shiny.R")
# source("data_visualization/Shiny_App/plot_continous_shiny.R")




source("dichotomous_functions.R")
source("continous_functions.R")

# Load plot generators
source("plot_dichotomous_shiny.R")
source("plot_continous_shiny.R")

# Model list - Dichotomous case
model_list_dic<-list("hill","gamma","logistic","log-probit","weibull",
                 "log-logistic","probit","multistage")

# Model list - Continous case
model_list_cont<-list("FUNL","hill","exp-3","exp-5","power")
fit_type<-list("mcmc", "mle", "laplace")




# Dichotomous case- Sample data input
mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 40,50,
                  50, 48,50),nrow=6,ncol=3,byrow=T)

# Example for continous case 
cont_dose<-c(0,0,0,0,0,18,18,18,18,20,20,20,20,30,30,30,30,35,35,35,35,39,39,39,39,39)
cont_response<- c(39.0,39,38.4,36.3,37.1,40.2,45.3,42.1,38.3,42.5,45.2,40.1,39.8,50.1,53.4,48.2,52.1,56.1,50.4,53.2,
                  55.2,55.1,59.1,56.3,52.9,53.7)

M<-matrix(nrow=length(cont_dose),ncol=2)
M[,1]<-cont_dose
M[,2]<-cont_response
M


ui<-navbarPage(title = "Toxic R - Interactive Plot V0.5", selected="Dichotomous Fitting",
               # Option 1. Select box option
               
               tabPanel("Dichotomous Fitting",  fluidPage(
                 sidebarLayout(
                   sidebarPanel(
                     
                     
                     conditionalPanel(condition="input.tabs =='Single Model'",
                                      helpText("Dichotomous Single Model"),
                                      
                                      selectInput(inputId="model", 
                                                  label= "Choose a model type",
                                                  choices=model_list_dic, 
                                                  selected = "gamma"),
                                      
                                      selectInput(inputId="fit_type", 
                                                  label= "Choose a fit type",
                                                  choices=fit_type, 
                                                  selected = "mcmc"),
                                      
                                      sliderInput(inputId="bmr_slide",
                                                  label="Choose a BMR level",
                                                  min=0,max=1,value=0.1)
                     ),
                     
                     conditionalPanel(condition="input.tabs =='Model Average'",
                                      helpText("Dichotomous Model Average"),
                                      
                                      # selectInput(inputId="model", 
                                      #             label= "Choose a model type",
                                      #             choices=model_list, 
                                      #             selected = "gamma"),
                                      # 
                                      selectInput(inputId="fit_type2",
                                                  label= "Choose a fit type",
                                                  choices=fit_type,
                                                  selected = "mcmc"),
                                      
                                      sliderInput(inputId="bmr_slide2",
                                                  label="Choose a BMR level",
                                                  min=0,max=1,value=0.1)
                     )
                     
                     
                   ),
                   mainPanel(
                     
                     tabsetPanel(id="tabs",
                                 tabPanel("Single Model",plotlyOutput(outputId = "dic_sing_plot")),
                                 tabPanel("Model Average",plotlyOutput(outputId = "dic_ma_plot"))
                                 
                                 
                     )
                     
                   )
                 )
                 
               )
               
               ),
               tabPanel("Continous Fitting",  fluidPage(
                 sidebarLayout(
                   sidebarPanel(
                     
                     
                     conditionalPanel(condition="input.tabs =='Single Model'",
                                      helpText("Continous Single Model"),
                                      
                                      selectInput(inputId="model3", 
                                                  label= "Choose a model type",
                                                  choices=model_list_cont, 
                                                  selected = "hill"),
                                      
                                      selectInput(inputId="fit_type3", 
                                                  label= "Choose a fit type",
                                                  choices=fit_type, 
                                                  selected = "mcmc"),
                                      
                                      sliderInput(inputId="bmr_slide3",
                                                  label="Choose a BMR level",
                                                  min=0,max=1,value=0.1)
                     ),
                     
                     conditionalPanel(condition="input.tabs =='Model Average'",
                                      helpText("Continous Model Average"),
                                      
                                       # selectInput(inputId="model2",
                                       #           label= "Choose a model type",
                                       #            choices=model_list_cont,
                                       #            selected = "mcmc"),
                                       # 
                                      selectInput(inputId="fit_type4",
                                                  label= "Choose a fit type",
                                                  choices=fit_type,
                                                  selected = "mcmc"),
                                      
                                      sliderInput(inputId="bmr_slide4",
                                                  label="Choose a BMR level",
                                                  min=0,max=1,value=0.1)
                     )
                     
                     
                   ),
                   mainPanel(
                     
                     tabsetPanel(id="tabs",
                                 tabPanel("Single Model",plotOutput(outputId = "cont_sing_plot")),
                                 tabPanel("Model Average",plotOutput(outputId = "cont_ma_plot"))
                                 
                                
                     )
                     
                   )
                 )
                 
               )
               
               )
)




server<- function (input,output){
  
  
  
  # Dichotomous Single Case
  output$dic_sing_plot<-renderPlotly({
    temp_fit = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = input$model,fit_type = input$fit_type, BMR = input$bmr_slide)
    
    #For MCMC  
    if (input$fit_type=="mcmc"){
      .plot.BMDdich_fit_MCMC(fit=temp_fit,fit_type=input$fit_type)
    }
    
    #For MCMC  
    else if (input$fit_type!="mcmc"){
      .plot.BMDdich_fit_maximized(fit=temp_fit,fit_type=input$fit_type)
    }
    
    
  })
  
  # Dichotomous Model Average Case
  output$dic_ma_plot<-renderPlotly({
    temp_fit = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = input$fit_type, BMR = input$bmr_slide2)
    .plot.BMDdichotomous_MA(A=temp_fit)
  })
  
  
  
  
  # Continous Single Case
  output$cont_sing_plot<-renderPlot({
    
    # Data input needs to be checke
    
    #temp_fit = single_continuous_fit(M[,1,drop=F],M[,2,drop=F],model_type = input$model2,fit_type = input$fit_type2, BMR = input$bmr_slide2)
    temp_fit2 <- single_continuous_fit(M[,1,drop=F],M[,2,drop=F],sstat = F, BMR = input$bmr_slide3 ,model_type=input$model3 ,distribution = "normal",fit_type = input$fit_type3)

    #For MCMC  
    if (input$fit_type2=="mcmc"){
      .plot.BMDcont_fit_MCMC(fit=temp_fit2,qprob=0.05)
    }
    
    #For MCMC  
    else if (input$fit_type2!="mcmc"){
      .plot.BMDcont_fit_maximized(fit=temp_fit2,qprob=0.05)
    }
    
  })  
  # Continous Model average case
  
  output$cont_ma_plot<-renderPlot({
    temp_fit2 = ma_continuous_fit(D=M[,1,drop=F],Y=M[,2,drop=F],model_list = NA, BMR = input$bmr_slide4 ,distribution = "normal",fit_type = input$fit_type4)
    .plot.BMDcontinuous_MA(temp_fit2,qprob=0.05)
  })  
  
  
}

shinyApp(ui=ui, server=server)

