library(shiny)
library(splines)
library(ggplot2)
library(scales)

# Define UI ----
ui <- fluidPage(
  titlePanel("Predicting adaptation and persistence with cubic splines"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      helpText("Fitting parameters:"),
      numericInput("df", "Degrees of freedom in spline:", 6),
      # numericInput("nx", "Phenotypic segments:", 1000),
      numericInput("bs", "Number of bootstraps:", 0),
      
      helpText("Population parameters:"),
      # numericInput("V_p", "Phenotypic variance:", 29.05, min = 0),
      numericInput("V_a", "Additive genetic variance in trait:", 2.62, min = 0),
      numericInput("rmax", "Maximum Malthusian growth rate:", 0.49, min = 0),
      numericInput("gT", "Generation time:", 1.81, min = 0),
      numericInput("B", "Environmental sensitivity (slope of optimum trait vs. environmental variable):", -5.30),
      numericInput("b", "Plasticity (slope of trait vs. environmental variable):", -4.98),
      
      helpText("Enter plotting ranges:"),
      numericInput("min_rate", "Minimum rate of change to plot:", -1),
      numericInput("max_rate", "Maximum rate of change to plot:", 1),
      numericInput("min_growth", "Minimum growth rate to plot:", -0.1),
      numericInput("max_growth", "Maximum growth rate to plot:", 0.5),
      
      helpText("Use test data (great tit laying data and spring temperature; Vedder et al 2013 PLoS Biol) or input your own."),  
      checkboxInput("test_data", "Use test data", TRUE),
      helpText("Input data [column 1 = phenotype, column 2 = fitness W=exp(r), e.g., expected number of offspring per generation]:"),  
      fileInput("file1", "Choose CSV file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      checkboxInput("header", "Does the CSV have a header?", TRUE)
      
    ),
    
    mainPanel(
      # tableOutput('contents'),
      plotOutput("fits"),
      plotOutput("lag"),
      plotOutput("crit")
    )
    
  )
)

# Define server logic ----
server <- function(input, output) {
      
  #function to resample data, with replacement, same number of data points as original
  resampler <- function(data) {
    n <- nrow(data)
    resample.rows <- sample(1:n,size=n,replace=TRUE)
    return(data[resample.rows,])
  }
  
  #function to fit gaussian for fitness as a function of phenotype (this is a quadratic in log fitness)
  quadratic.estimator <- function(data, x) {
    quad_fit <- nls(fitness ~ a * exp( - (b - phenotype)^2 / (2*c^2)), data, start = list(a = max(data$fitness), b = data$phenotype[which.max(data$fitness)], c = 10))
    return(predict(quad_fit, list(phenotype = x))) #predicted values over x values
  }
  
  #function to get quadratic confidence intervals
  quad.cis <- function(data, x, B, alpha = 0.05) {
    quad.main <- quadratic.estimator(data, x) #fit on original data
    if(B == 0) {
      #if no bootstraps then just make the confidence interval of zero width
      cis.lower <- quad.main #lower CI
      cis.upper <- quad.main #upper CI      
    } else {
      quad.boots <- replicate(B, quadratic.estimator(resampler(data), x)) #resample B times
      cis.lower <- 2*quad.main - apply(quad.boots,1,quantile,probs=1-alpha/2) #lower CI
      cis.lower[which(cis.lower<0)] <- 0 #make all negative values zeros
      cis.upper <- 2*quad.main - apply(quad.boots,1,quantile,probs=alpha/2) #upper CI
    }
    return(list(main.curve = quad.main, lower.ci = cis.lower, upper.ci = cis.upper, x = x))
  }
  
  #function to fit spline
  spline.estimator <- function(data, x) {
    ss_fit <- smooth.spline(x = data$phenotype, y = data$fitness, df = input$df)
    pred <- predict(ss_fit, x = x)$y #predict over x
    pred[which(pred<0)] <- 0 #make all negative values zeros
    return(pred) # non-negative predicted values
  }
  
  #function to get spline confidence intervals
  spline.cis <- function(data,x,B,alpha=0.05) {
    spline.main <- spline.estimator(data,x) #fit on original data
    if(B == 0) {
      #if no bootstraps then just make the confidence interval of zero width
      cis.lower <- spline.main #lower CI
      cis.upper <- spline.main #upper CI         
    } else {
      spline.boots <- replicate(B,spline.estimator(resampler(data),x)) #resample B times
      cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2) #lower CI
      cis.lower[which(cis.lower<0)] <- 0 #make all negative values zeros
      cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2) #upper CI      
    }
    return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,x=x))
  }  

  #function to get derivative of quadratic
  quad.deriv <- function(data,x) {
    quad.main <- quadratic.estimator(data,x)
    quad.log <- log(quad.main)
    return(diff(quad.log)/diff(x))
  }
  
  #function to get derivative of spline
  spline.deriv <- function(data,x) {
    spline.main <- spline.estimator(data,x)
    spline.log <- log(spline.main)
    return(diff(spline.log)/diff(x))
  }
  
  #function to get quadratic derivative confidence intervals
  quad.deriv.cis <- function(data,x,B,alpha=0.05) {
    quad.main <- quad.deriv(data,x)
    if(B == 0) {
      cis.lower <- quad.main
      cis.upper <- quad.main
    } else {
      quad.boots <- replicate(B,quad.deriv(resampler(data),x))
      cis.lower <- 2*quad.main - apply(quad.boots,1,quantile,probs=1-alpha/2,na.rm=TRUE)
      cis.upper <- 2*quad.main - apply(quad.boots,1,quantile,probs=alpha/2,na.rm=TRUE)
    }
    dX <- rowMeans(embed(x,2)) #centred x values      
    return(list(main.curve=quad.main,lower.ci=cis.lower,upper.ci=cis.upper,
                x=dX))
  }
  
  #function to get spline derivative confidence intervals
  spline.deriv.cis <- function(data,x,B,alpha=0.05) {
    spline.main <- spline.deriv(data,x)
    if(B ==0) {
      cis.lower <- spline.main
      cis.upper <- spline.main     
    } else {
      spline.boots <- replicate(B,spline.deriv(resampler(data),x))
      cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2,na.rm=TRUE)
      cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2,na.rm=TRUE)      
    }
    dX <- rowMeans(embed(x,2)) #centred x values
    return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
                x=dX))
  }
    
  params <- reactiveValues()
  
  observe({
    
    #degreee to fit
    params$df <- input$df #degrees of freedom in spline fit
    # params$nx <- input$nx #number of phenotypic segments for derivatives and plots
    params$nx <- 1000
    params$bs <- input$bs
    
    #parameters
    # params$V_p <- input$V_p #phenotypic variance
    params$V_a <- input$V_a #additive genetic variance
    params$rmax <- input$rmax #max malthusian growth rate
    params$gT <- input$gT #generation time (avg age of mothers)
    params$B <- input$B #environmental sensitivity (slope of optimum phenotype vs environmental variable)
    params$b <- input$b #plasticity (slope of expressed phenotype vs environmental variable)
    
    #plot limits
    params$min_rate = input$min_rate
    params$max_rate = input$max_rate
    params$min_growth = input$min_growth
    params$max_growth = input$max_growth
    
  })
  
  myData <- reactive({
    
    if(input$test_data) {
      #if use test data then use phenotype and fitness data pulled from Vedder et al. 2013 PLoS Biol. 11(7):1001605, supplementary figure 1
      phenotype <- c(-18, -14, -10, -6, -2, 2, 6, 10, 14, 18) #phenotype (mean of bin)
      fitness <- c(0.18, 0.95, 0.8, 0.78, 0.76, 0.65, 0.45, 0.4, 0.35, 0.2) #expected fitness of phenotype (mean fitness within bin)
      data <- data.frame(phenotype, fitness) #combine to data frame
    } else {
      # req(input$file1)
      inFile <- input$file1
      data <- read.csv(inFile$datapath, header=input$header)
    }

    data[,3] <- log(data[,2]) #add log fitness as column
    names(data) <- c('phenotype', 'fitness', 'log_fitness')
    data
    
  })

  # output$contents <- renderTable({
  #   myData()
  # })  
    
  #enter phenotype and fitness data
  # phenotype <- c(-18, -14, -10, -6, -2, 2, 6, 10, 14, 18) #phenotype (mean of bin)
  # # if(input$flip){phenotype <- phenotype * -1} #flip sign of phenotype if interested in decreasing environment
  # fitness <- c(0.18, 0.95, 0.8, 0.78, 0.76, 0.65, 0.45, 0.4, 0.35, 0.2) #expected fitness of phenotype (mean fitness within bin)
  # data <- data.frame(phenotype, fitness) #combine to data frame
  # data$log_fitness <- log(data$fitness) #add log fitness as column
  # data <- read.csv(file="vedder_data.csv", header=TRUE, sep=",")
  # data$log_fitness <- log(data$fitness) #add log fitness as column

  output$fits <- renderPlot({
    
    data <- myData()
    
    # phenotype <- data
    # #x-values for plotting and derivatives
    x <- with(data, seq(min(phenotype), max(phenotype), (max(phenotype)-min(phenotype))/params$nx)) #divide up the x-axis

    # #fit a quadratic and get strength of selection and optimum
    # nlsfit <- nls(log_fitness ~ a - (b - phenotype)^2 / (2*c^2), data, start=list(a=0,b=1,c=0.1)) #fit quadratic over data to get QG parameters
    # quad_fit <- predict(nlsfit, list(phenotype=x)) #fit over chosen phenotypic segments
    # # gamma <- (input$V_p + as.numeric(coef(nlsfit)[3])^2)^(-1) #effective strength of selection (gamma)
    # # gamma <- (input$V_p + 11.62^2)^(-1) #effective strength of selection in vedder (gamma)
    # opt_quad <- as.numeric(coef(nlsfit)[2]) #phenotype that maximizes fitness in quadratic
    # 
    # #fit a cubic spline and get optimum
    # ss <- smooth.spline(data$phenotype, data$log_fitness, df=params$df)
    # ss_fit <- predict(ss, x) #fit over chosen phenotypic segments
    # opt_spline <- x[which.max(ss_fit$y)] #estimate phenotype that maximizes fitness
    # 
    # #examine fits visually
    # plot_data <- data.frame(x, quad_fit, ss_fit$y) #combine fits and phenotype
    # ggplot() +
    #   geom_point(data = data, aes(x = phenotype, y = log_fitness), colour="black") +
    #   geom_line(data = plot_data, aes(x = x, y = quad_fit), colour="red") +
    #   geom_line(data = plot_data, aes(x = x, y = ss_fit.y), colour="blue") +
    #   labs(x = 'phenotype', y = 'log fitness') +
    #   theme_bw() +
    #   theme(panel.border= element_blank()) + 
    #   annotate("text", x = max(x), y = max(data$log_fitness), label = "quadratic", hjust=1, vjust = 0, color="red") +
    #   annotate("text", x = max(x), y = max(data$log_fitness), label = "spline", hjust=1, vjust = 2, color="blue")
  
    #run bootstrap
    q.cis <- as.data.frame(quad.cis(data, x, B=params$bs, alpha=0.05))
    sp.cis <- as.data.frame(spline.cis(data, x, B=params$bs, alpha=0.05))
    
    #plot on log scale with confidence intervals
    ggplot() +
      geom_point(data = data, aes(x = phenotype, y = fitness), colour = "gray") +
      geom_line(data = q.cis, aes(x = x, y = main.curve), colour = "red", size = 1) +
      # geom_line(data = q.cis, aes(x = x, y = lower.ci), colour="red", lty=2) +
      # geom_line(data = q.cis, aes(x = x, y = upper.ci), colour="red", lty=2) +
      geom_ribbon(data = q.cis, aes(x = x, ymin = lower.ci, ymax = upper.ci), fill = "red", alpha = "0.5") +
      geom_line(data = sp.cis, aes(x = x, y = main.curve), colour="blue", size=1) +
      # geom_line(data = sp.cis, aes(x = x, y = lower.ci), colour="blue", lty=2) +
      # geom_line(data = sp.cis, aes(x = x, y = upper.ci), colour="blue", lty=2) +
      geom_ribbon(data = sp.cis, aes(x = x, ymin = lower.ci, ymax = upper.ci), fill = "blue", alpha = "0.5") +
      scale_y_continuous(trans = log_trans(), breaks = c(1, 10, 100, 1000, 10000)) +
      labs(x = 'phenotype', y = 'fitness (log scale)') + 
      theme_bw() + 
      theme(panel.border = element_blank()) + 
      annotate("text", x = max(x), y = max(data$fitness), label = "quadratic", hjust = 1, vjust = 0, color = "red") +
      annotate("text", x = max(x), y = max(data$fitness), label = "spline", hjust = 1, vjust = 2, color = "blue")
      
  })
   
  output$lag <- renderPlot({
    
    data <- myData()
    x <- with(data, seq(min(phenotype), max(phenotype), (max(phenotype)-min(phenotype))/params$nx)) #divide up the x-axis
    
    #constant to get rate of evo from selction gradient
    evo_scale = params$V_a / ((params$B - params$b) * params$gT) 
    
    # phenotype <- data[,1]
    # log_fitness <- data[,2]
    # 
    # #x-values for plotting and derivatives
    # x <- seq(min(phenotype), max(phenotype), (max(phenotype)-min(phenotype))/params$nx) #divide up the x-axis
    # 
    # #fit a quadratic and get strength of selection and optimum
    # nlsfit <- nls(log_fitness ~ a - (b - phenotype)^2 / (2*c^2), data, start=list(a=0,b=1,c=0.1)) #fit quadratic over data to get QG parameters
    # quad_fit <- predict(nlsfit, list(phenotype=x)) #fit over chosen phenotypic segments
    # # gamma <- (input$V_p + as.numeric(coef(nlsfit)[3])^2)^(-1) #effective strength of selection (gamma)
    # # gamma <- (input$V_p + 11.62^2)^(-1) #effective strength of selection in vedder (gamma)
    # opt_quad <- as.numeric(coef(nlsfit)[2]) #phenotype that maximizes fitness in quadratic
    # 
    # #fit a cubic spline and get optimum
    # ss <- smooth.spline(data$phenotype, data$log_fitness, df=params$df)
    # ss_fit <- predict(ss, x) #fit over chosen phenotypic segments
    # opt_spline <- x[which.max(ss_fit$y)] #estimate phenotype that maximizes fitness
    # 
    # #take derivatives
    # dY_quad <- diff(quad_fit)/diff(x) #slope of quadratic
    # dY_ss <- diff(ss_fit$y)/diff(x) #slope of spline
    # dX <- rowMeans(embed(x,2)) #centred x values
    # # dY_quad <- 1/(as.numeric(coef(nlsfit)[3])^2 + input$V_p)*(opt_quad - dX) #selection gradient as determined by theory
    # 
    # #maxima in derivative over positive lags
    # pos_lags <- which(opt_spline - dX > 0) #index for positive lags
    # pos_dY_ss <- dY_ss[pos_lags] #derivatives over positive lags
    # pos_maxs <- which(diff(sign(diff(pos_dY_ss)))==-2) + 1 #local maxima in derivative over positive lags
    # pos_nmaxs <- length(pos_maxs) #number of maxima
    # 
    # #minima in derivative over negative lags
    # neg_lags <- which(opt_spline - dX < 0) #index for negative lags
    # neg_dY_ss <- dY_ss[neg_lags] #derivatives over negative lags
    # neg_mins <- which(diff(sign(diff(neg_dY_ss)))==2) + 1 #local minima in derivative over negative lags
    # neg_nmins <- length(neg_mins) #number of minima
    # 
    # #constant to get rate of evo from selction gradient
    # evo_scale = params$V_a / ((params$B - params$b) * params$gT) 
    # 
    # #steady state lag as a function of the rate of environmental change
    # plot = 
    #   ggplot() +
    #   geom_line(aes(x = evo_scale * dY_quad, y = opt_quad - dX), colour="red") +
    #   labs(x = 'rate of environmental change', y = 'steady-state phenotypic lag') +
    #   theme_bw() + 
    #   theme(panel.border= element_blank()) +
    #   geom_hline(yintercept = 0) + 
    #   geom_vline(xintercept = 0) + 
    #   annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(opt_quad - dX, opt_spline - dX)), label = "quadratic", hjust=1, vjust = 0, color="red") +
    #   annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(opt_quad - dX, opt_spline - dX)), label = "spline", hjust=1, vjust = 2, color="blue")
    # 
    # #plot spline over positive lags
    # if(pos_nmaxs == 0) {
    #   #if no maxs then just plot solid curve over all lags
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * pos_dY_ss, y = opt_spline - dX[pos_lags]), colour="blue")
    # } else {
    #   #if some maxs first plot dashed curve over all lags and solid curve from lag of zero to first max
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * pos_dY_ss, y = opt_spline - dX[pos_lags]), colour="blue", lty=2) + 
    #     geom_path(aes(x = evo_scale * pos_dY_ss[length(dY_ss[pos_lags]):pos_maxs[pos_nmaxs]], y = opt_spline - dX[pos_lags][length(dY_ss[pos_lags]):pos_maxs[pos_nmaxs]]), colour="blue") 
    #   keep_going <- TRUE
    #   current_index <- pos_maxs[pos_nmaxs]
    #   while(keep_going) {
    #     if(any(pos_dY_ss[current_index] < pos_dY_ss & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags])) {
    #       #if there is a greater lag that gives a greater derivative, set this to be the next stable lag
    #       next_index <- max(which(pos_dY_ss[current_index] < dY_ss[pos_lags] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags]))
    #       if(any(pos_dY_ss[current_index] < pos_dY_ss[pos_maxs] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags][pos_maxs])) {
    #         #if there is another max at a greater lag, plot the solid line only up to there
    #         last_index <- pos_maxs[max(which(pos_dY_ss[current_index] < pos_dY_ss[pos_maxs] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags][pos_maxs]))]
    #         current_index <- last_index
    #       } else {
    #         #if there are no more maxs then plot all the way to the end of the data and stop the loop
    #         last_index <- 1
    #         keep_going <- FALSE
    #       }
    #       #make temporary data frame with rate of enviro change and steady-state lag over range of interest (this lets ggplot know we want to plot a new layer each time)
    #       gg.data <- data.frame(x=evo_scale * pos_dY_ss[next_index:last_index], y=opt_spline - dX[pos_lags][next_index:last_index])
    #       plot = plot +
    #         geom_path(data = gg.data, aes(x = x, y = y), colour="blue")
    #     } else {
    #       #if there are no more maxs then stop the loop
    #       keep_going <- FALSE
    #     }
    #   }  
    # } 
    # #plot spline over negative lags
    # if(neg_nmins == 0) {
    #   #if no mins then just plot solid curve over all lags
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * neg_dY_ss, y = opt_spline - dX[neg_lags]), colour="blue")
    # } else {
    #   #if some mins first plot dashed curve over all lags and solid curve to first min
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * neg_dY_ss, y = opt_spline - dX[neg_lags]), colour="blue", lty=2) + 
    #     geom_path(aes(x = evo_scale * neg_dY_ss[1:neg_mins[1]], y = opt_spline - dX[neg_lags][1:neg_mins[1]]), colour="blue") 
    #   keep_going <- TRUE
    #   current_index <- neg_mins[1]
    #   while(keep_going) {
    #     if(any(neg_dY_ss[current_index] > neg_dY_ss & opt_spline - dX[neg_lags][current_index] > opt_spline - dX[neg_lags])) {
    #       next_index <- min(which(neg_dY_ss[current_index] > neg_dY_ss & opt_spline - dX[neg_lags][current_index] > opt_spline - dX[neg_lags]))
    #       if(any(neg_dY_ss[next_index] > neg_dY_ss[neg_mins] & opt_spline - dX[neg_lags][next_index] > opt_spline - dX[neg_lags][neg_mins])) {
    #         last_index <- neg_mins[min(which(neg_dY_ss[next_index] > neg_dY_ss[neg_mins] & opt_spline - dX[neg_lags][next_index] > opt_spline - dX[neg_lags][neg_mins]))]
    #         current_index <- last_index
    #       } else {
    #         last_index <- length(neg_lags)
    #         keep_going <- FALSE
    #       }
    #       gg.data <- data.frame(x=evo_scale * neg_dY_ss[next_index:last_index], y=opt_spline - dX[neg_lags][next_index:last_index])
    #       plot = plot +
    #         geom_path(data = gg.data, aes(x = x, y = y), colour="blue")
    #     } else {
    #       keep_going <- FALSE
    #     }
    #   }  
    # }
    # plot
    
    #perform bootstrap
    q.d.cis <- as.data.frame(quad.deriv.cis(data, x, B=params$bs, alpha=0.05))
    sp.d.cis <- as.data.frame(spline.deriv.cis(data, x, B=params$bs, alpha=0.05))
    
    #get optimal trait values
    quad_opt <- x[which.max(quadratic.estimator(data,x))]
    spline_opt <- x[which.max(spline.estimator(data,x))]
    
    #maxima in derivative over positive lags
    pos_lags <- which(spline_opt - sp.d.cis$x > 0) #index for positive lags
    pos_dY_ss <- sp.d.cis$main.curve[pos_lags] #derivatives over positive lags
    pos_maxs <- which(diff(sign(diff(pos_dY_ss)))==-2) + 1 #local maxima in derivative over positive lags
    pos_nmaxs <- length(pos_maxs) #number of maxima
    
    #minima in derivative over negative lags
    neg_lags <- which(spline_opt - sp.d.cis$x < 0) #index for negative lags
    neg_dY_ss <- sp.d.cis$main.curve[neg_lags] #derivatives over negative lags
    neg_mins <- which(diff(sign(diff(neg_dY_ss)))==2) + 1 #local minima in derivative over negative lags
    neg_nmins <- length(neg_mins) #number of minima
    
    #steady state lag as a function of the rate of environmental change
    xmin = min(c(quad_opt - q.d.cis$x, spline_opt - sp.d.cis$x))
    xmax = max(c(quad_opt - q.d.cis$x, spline_opt - sp.d.cis$x))
    ymin = params$min_rate
    ymax = params$max_rate
    plot = 
      ggplot() +
      geom_line(data = q.d.cis, aes(x = quad_opt - x, y = evo_scale * main.curve), colour="red", size=1) +
      geom_ribbon(data = q.d.cis, aes(x = quad_opt - x, ymin = evo_scale * lower.ci, ymax = evo_scale * upper.ci), fill="red", alpha="0.5") +
      # geom_line(data = sp.d.cis, aes(x = spline_opt - x, y = evo_scale * main.curve), colour="blue", size=1) +
      geom_ribbon(data = sp.d.cis, aes(x = spline_opt - x, ymin = evo_scale * lower.ci, ymax = evo_scale * upper.ci), fill="blue", alpha="0.5") +
      labs(x = 'steady-state phenotypic lag', y = 'rate of environmental change') +
      theme_bw() + 
      theme(panel.border= element_blank()) +
      geom_hline(yintercept = 0) + 
      geom_vline(xintercept = 0) + 
      ylim(ymin,ymax) +
      annotate("text", x = xmax, y = ymax, label = "quadratic", hjust=1, vjust = 0, color="red") +
      annotate("text", x = xmax, y = ymax, label = "spline", hjust=1, vjust = 2, color="blue") +
      coord_flip()
    #plot spline over positive lags
    if(pos_nmaxs == 0) {
      #if no maxs then just plot solid curve over all lags
      plot = plot + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[pos_lags], y = evo_scale * pos_dY_ss), colour="blue", size=1)
    } else {
      #if some maxs first plot dashed curve over all lags and solid curve from lag of zero to first max
      plot = plot + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[pos_lags], y = evo_scale * pos_dY_ss), colour="blue", lty=2) + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[pos_lags][length(pos_lags):pos_maxs[pos_nmaxs]], y = evo_scale * pos_dY_ss[length(pos_lags):pos_maxs[pos_nmaxs]]), colour="blue", size=1) 
      keep_going <- TRUE
      current_index <- pos_maxs[pos_nmaxs]
      while(keep_going) {
        if(any(pos_dY_ss[current_index] < pos_dY_ss & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags])) {
          #if there is a greater lag that gives a greater derivative, set this to be the next stable lag
          potentials <- which(pos_dY_ss[current_index] < pos_dY_ss & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags])
          next_index <- potentials[which.min(pos_dY_ss[potentials])]
          if(any(pos_dY_ss[current_index] < pos_dY_ss[pos_maxs] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags][pos_maxs])) {
            #if there is another max at a greater lag, plot the solid line only up to there
            last_index <- pos_maxs[max(which(pos_dY_ss[current_index] < pos_dY_ss[pos_maxs] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags][pos_maxs]))]
            current_index <- last_index
          } else {
            #if there are no more maxs then plot all the way to the end of the data and stop the loop
            last_index <- 1
            keep_going <- FALSE
          }
          #make temporary data frame with rate of enviro change and steady-state lag over range of interest (this lets ggplot know we want to plot a new layer each time)
          gg.data <- data.frame(x = spline_opt - sp.d.cis$x[pos_lags][next_index:last_index], y = evo_scale * pos_dY_ss[next_index:last_index])
          plot = plot +
            geom_path(data = gg.data, aes(x = x, y = y), colour="blue", size=1)
        } else {
          #if there are no more maxs then stop the loop
          keep_going <- FALSE
        }
      }  
    }
    #plot spline over negative lags
    if(neg_nmins == 0) {
      #if no mins then just plot solid curve over all lags
      plot = plot + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[neg_lags], y = evo_scale * neg_dY_ss), colour="blue", size=1)
    } else {
      #if some mins first plot dashed curve over all lags and solid curve from lag of zero to first min
      plot = plot + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[neg_lags], y = evo_scale * neg_dY_ss), colour="blue", lty=2) + 
        geom_path(aes(x = spline_opt - sp.d.cis$x[neg_lags][1:neg_mins[1]], y = evo_scale * neg_dY_ss[1:neg_mins[1]]), colour="blue", size=1) 
      keep_going <- TRUE
      current_index <- neg_mins[1]
      while(keep_going) {
        if(any(neg_dY_ss[current_index] > neg_dY_ss & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags])) {
          #if there is a greater lag that gives a greater derivative, set this to be the next stable lag
          potentials <- which(neg_dY_ss[current_index] > neg_dY_ss & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags])
          next_index <- potentials[which.max(neg_dY_ss[potentials])]
          if(any(neg_dY_ss[current_index] > neg_dY_ss[neg_mins] & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags][neg_mins])) {
            #if there is another max at a greater lag, plot the solid line only up to there
            last_index <- neg_mins[max(which(neg_dY_ss[current_index] > neg_dY_ss[neg_mins] & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags][neg_mins]))]
            current_index <- last_index
          } else {
            #if there are no more maxs then plot all the way to the end of the data and stop the loop
            last_index <- length(neg_lags)
            keep_going <- FALSE
          }
          #make temporary data frame with rate of enviro change and steady-state lag over range of interest (this lets ggplot know we want to plot a new layer each time)
          gg.data <- data.frame(x = spline_opt - sp.d.cis$x[neg_lags][next_index:last_index], y = evo_scale * neg_dY_ss[next_index:last_index])
          plot = plot +
            geom_path(data = gg.data, aes(x = x, y = y), colour="blue", size=1)
        } else {
          #if there are no more maxs then stop the loop
          keep_going <- FALSE
        }
      }  
    } 
    plot
    
  })

  output$crit <- renderPlot({
    
    data <- myData()
    # phenotype <- data[,1]
    # log_fitness <- data[,2]
    # 
    # #x-values for plotting and derivatives
    # x <- seq(min(phenotype), max(phenotype), (max(phenotype)-min(phenotype))/params$nx) #divide up the x-axis
    # 
    # #fit a quadratic and get strength of selection and optimum
    # nlsfit <- nls(log_fitness ~ a - (b - phenotype)^2 / (2*c^2), data, start=list(a=0,b=1,c=0.1)) #fit quadratic over data to get QG parameters
    # quad_fit <- predict(nlsfit, list(phenotype=x)) #fit over chosen phenotypic segments
    # # gamma <- (input$V_p + as.numeric(coef(nlsfit)[3])^2)^(-1) #effective strength of selection (gamma)
    # # gamma <- (input$V_p + 11.62^2)^(-1) #effective strength of selection in vedder (gamma)
    # opt_quad <- as.numeric(coef(nlsfit)[2]) #phenotype that maximizes fitness in quadratic
    # 
    # #fit a cubic spline and get optimum
    # ss <- smooth.spline(data$phenotype, data$log_fitness, df=params$df)
    # ss_fit <- predict(ss, x) #fit over chosen phenotypic segments
    # opt_spline <- x[which.max(ss_fit$y)] #estimate phenotype that maximizes fitness
    # 
    # #take derivatives
    # dY_quad <- diff(quad_fit)/diff(x) #slope of quadratic
    # dY_ss <- diff(ss_fit$y)/diff(x) #slope of spline
    # dX <- rowMeans(embed(x,2)) #centred x values
    # # dY_quad <- 1/(as.numeric(coef(nlsfit)[3])^2 + input$V_p)*(opt_quad - dX) #selection gradient as determined by theory
    # 
    # #maxima in derivative over positive lags
    # pos_lags <- which(opt_spline - dX > 0) #index for positive lags
    # pos_dY_ss <- dY_ss[pos_lags] #derivatives over positive lags
    # pos_maxs <- which(diff(sign(diff(pos_dY_ss)))==-2) + 1 #local maxima in derivative over positive lags
    # pos_nmaxs <- length(pos_maxs) #number of maxima
    # 
    # #minima in derivative over negative lags
    # neg_lags <- which(opt_spline - dX < 0) #index for negative lags
    # neg_dY_ss <- dY_ss[neg_lags] #derivatives over negative lags
    # neg_mins <- which(diff(sign(diff(neg_dY_ss)))==2) + 1 #local minima in derivative over negative lags
    # neg_nmins <- length(neg_mins) #number of minima
    # 
    # #constant to get rate of evo from selction gradient
    # evo_scale = params$V_a / ((params$B - params$b) * params$gT) 
    # 
    # #calculate the critical rate (or minimum positive growth rate) for positive rates of environmental change
    # pred_quad <- predict(nlsfit, list(phenotype=dX)) #list of quad predictions of fitness over phenotype that we use to calculate selection gradient
    # quad_max = max(pred_quad) #fitness of optimum in quadratic
    # pos_grow_quad <- which(params$rmax - (quad_max - pred_quad) / params$gT > 0) #list of rows that have positive growth
    # crit_rate <- max(evo_scale * dY_quad[pos_grow_quad]) #potential critical rate (max rate of evo with positive growth)
    # crit_index <- which.max(evo_scale * dY_quad[pos_grow_quad]) #get the index of the potential critical rate
    # crit_growth <- params$rmax - (quad_max - pred_quad[pos_grow_quad][crit_index]) / params$gT #get the growth rate at the potential critical rate
    # neg_rates_quad <- any(evo_scale * dY_quad > 0 & params$rmax - (quad_max - pred_quad) / params$gT < 0) #are there positive rates of change with negative growth rates (ie can we esimtate a critical rate or does the data not extend far enough)
    # 
    # #calculate the critical rate (or minimum positive growth rate) for negative rates of environmental change
    # crit_rate_min <- min(evo_scale * dY_quad[pos_grow_quad]) #potential critical rate (min rate of evo with positive growth)
    # crit_index_min <- which.min(evo_scale * dY_quad[pos_grow_quad]) #get the index of the potential critical rate
    # crit_growth_min <- params$rmax - (quad_max - pred_quad[pos_grow_quad][crit_index_min]) / params$gT #get the growth rate at the potential critical rate
    # neg_rates_quad_min <- any(evo_scale * dY_quad < 0 & params$rmax - (quad_max - pred_quad) / params$gT < 0) #are there negative rates of evo with negative growth rates (ie can we esimtate a critical rate or does the data not extend far enough)
    # 
    # #calculate the tipping point (or minimum positive growth rate) for positive rates of environmental change
    # pred_ss <- predict(ss, dX)$y #list of spline predictions of fitness over phenotype that we use to calculate selection gradient
    # ss_max = max(pred_ss) #fitness of optimum in spline
    # pos_grow <- which(params$rmax - (ss_max - pred_ss) / params$gT > 0) #list of rows that have positive growth
    # tp_pt_rate <- max(evo_scale * dY_ss[pos_grow]) #potential tipping point
    # tp_index <- which.max(evo_scale * dY_ss[pos_grow]) #get the index
    # tp_pt_growth <- params$rmax - (ss_max - pred_ss[pos_grow][tp_index]) / params$gT #growth rate at tipping point
    # neg_rates_ss <- any(evo_scale * dY_ss > 0 & params$rmax - (ss_max - pred_ss) / params$gT < 0) #are there positive rates of change with negative growth rate (ie can we esimtate a tipping point)
    # 
    # #calculate the tipping point (or minimum positive growth rate) for negative rates of environmental change
    # tp_pt_rate_min <- min(evo_scale * dY_ss[pos_grow]) #potential tipping point
    # tp_index_min <- which.min(evo_scale * dY_ss[pos_grow]) #get the index
    # tp_pt_growth_min <- params$rmax - (ss_max - pred_ss[pos_grow][tp_index_min]) / params$gT #growth rate at tipping point
    # neg_rates_ss_min <- any(evo_scale * dY_ss < 0 & params$rmax - (ss_max - pred_ss) / params$gT < 0) #are there negative rates of evo with negative growth rates (ie can we esimtate a tipping point)
    # 
    # #rate of popn growth as a function of the rate of environmental change
    # plot = 
    #   ggplot() +
    #   geom_line(aes(x = evo_scale * dY_quad, y = input$rmax - (quad_max - pred_quad)/input$gT), colour="red") +
    #   labs(x = 'rate of environmental change', y = 'long-term popn malthusian growth rate') +
    #   # xlim(input$min_rate, input$max_rate) +
    #   # ylim(input$min_growth, input$max_growth) +
    #   annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "~underline('max sustainable rates')", parse = TRUE, hjust=1, vjust = 0, color="black") +
    #   annotate("text", x = min(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "~underline('min sustainable rates')", parse = TRUE, hjust=0, vjust = 0, color="black") +
    #   theme_bw() + 
    #   theme(panel.border= element_blank()) +
    #   geom_hline(yintercept = 0) + 
    #   geom_vline(xintercept = 0)
    # #plot spline over positive lags
    # if(pos_nmaxs == 0) {
    #   #if no maxs then just plot solid curve over all lags
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * dY_ss[pos_lags], y = input$rmax - (ss_max - pred_ss[pos_lags])/input$gT), colour="blue")
    # } else {
    #   #if some maxs first plot dashed curve over all lags and solid curve from lag of zero to first max
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * dY_ss[pos_lags], y = input$rmax - (ss_max - pred_ss[pos_lags])/input$gT), colour="blue", lty=2) + 
    #     geom_path(aes(x = evo_scale * dY_ss[pos_lags][length(dY_ss[pos_lags]):pos_maxs[pos_nmaxs]], y = input$rmax - (ss_max - pred_ss[pos_lags][length(dY_ss[pos_lags]):pos_maxs[pos_nmaxs]])/input$gT), colour="blue") 
    #   keep_going <- TRUE
    #   current_index <- pos_maxs[pos_nmaxs]
    #   while(keep_going) {
    #     if(any(dY_ss[pos_lags][current_index] < dY_ss[pos_lags] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags])) {
    #       #if there is a greater lag that gives a greater derivative, set this to be the next stable lag
    #       next_index <- max(which(dY_ss[pos_lags][current_index] < dY_ss[pos_lags] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags]))
    #       if(any(dY_ss[pos_lags][current_index] < dY_ss[pos_lags][pos_maxs] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags][pos_maxs])) {
    #         #if there is another max at a greater lag, plot the solid line only up to there
    #         last_index <- pos_maxs[max(which(dY_ss[pos_lags][current_index] < dY_ss[pos_lags][pos_maxs] & opt_spline - dX[pos_lags][current_index] < opt_spline - dX[pos_lags][pos_maxs]))]
    #         current_index <- last_index
    #       } else {
    #         #if there are no more maxs then plot all the way to the end of the data and stop the loop
    #         last_index <- 1
    #         keep_going <- FALSE
    #       }
    #       #make temporary data frame with rate of enviro change and steady-state lag over range of interest (this lets ggplot know we want to plot a new layer each time)
    #       gg.data <- data.frame(x=evo_scale * dY_ss[pos_lags][next_index:last_index], y=input$rmax - (ss_max - pred_ss[pos_lags][next_index:last_index])/input$gT)
    #       plot = plot +
    #         geom_path(data = gg.data, aes(x = x, y = y), colour="blue")
    #     } else {
    #       #if there are no more maxs then stop the loop
    #       keep_going <- FALSE
    #     }
    #   }  
    # } 
    # #plot spline over negative lags
    # if(neg_nmins == 0) {
    #   #if no mins then just plot solid curve over all lags
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * dY_ss[neg_lags], y = input$rmax - (ss_max - pred_ss[neg_lags])/input$gT), colour="blue")
    # } else {
    #   #if some mins first plot dashed curve over all lags and solid curve to first min
    #   plot = plot + 
    #     geom_path(aes(x = evo_scale * dY_ss[neg_lags], y = input$rmax - (ss_max - pred_ss[neg_lags])/input$gT), colour="blue", lty=2) + 
    #     geom_path(aes(x = evo_scale * dY_ss[neg_lags][1:neg_mins[1]], y = input$rmax - (ss_max - pred_ss[neg_lags][1:neg_mins[1]])/input$gT), colour="blue") 
    #   keep_going <- TRUE
    #   current_index <- neg_mins[1]
    #   while(keep_going) {
    #     if(any(dY_ss[neg_lags][current_index] > dY_ss[neg_lags] & opt_spline - dX[neg_lags][current_index] > opt_spline - dX[neg_lags])) {
    #       next_index <- min(which(dY_ss[neg_lags][current_index] > dY_ss[neg_lags] & opt_spline - dX[neg_lags][current_index] > opt_spline - dX[neg_lags]))
    #       if(any(dY_ss[neg_lags][next_index] > dY_ss[neg_lags][neg_mins] & opt_spline - dX[neg_lags][next_index] > opt_spline - dX[neg_lags][neg_mins])) {
    #         last_index <- neg_mins[min(which(dY_ss[neg_lags][next_index] > dY_ss[neg_lags][neg_mins] & opt_spline - dX[neg_lags][next_index] > opt_spline - dX[neg_lags][neg_mins]))]
    #         current_index <- last_index
    #       } else {
    #         last_index <- length(neg_lags)
    #         keep_going <- FALSE
    #       }
    #       gg.data <- data.frame(x=evo_scale * dY_ss[neg_lags][next_index:last_index], y=input$rmax - (ss_max - pred_ss[neg_lags][next_index:last_index])/input$gT)
    #       plot = plot +
    #         geom_path(data = gg.data, aes(x = x, y = y), colour="blue")
    #     } else {
    #       keep_going <- FALSE
    #     }
    #   }  
    # }
    # #print and plot critical rate for positive change if possible
    # if(neg_rates_quad) { #if we can estimate a critical rate
    #   plot = plot + 
    #     annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = paste("quadratic:", round(crit_rate,3)), hjust=1, vjust = 2, color="red") +
    #     geom_point(aes(x=crit_rate, y=0), color="red", size = 3) 
    # } else {
    #   plot = plot + 
    #     annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "quadratic: NA", hjust=1, vjust = 2, color="red")
    # }
    # #print and plot critical rate for negative change if possible
    # if(neg_rates_quad_min) { #if we can estimate a critical rate
    #   plot = plot + 
    #     annotate("text", x = min(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = paste("quadratic:", round(crit_rate_min,3)), hjust=0, vjust = 2, color="red") +
    #     geom_point(aes(x=crit_rate_min, y=0), color="red", size = 3) 
    # } else {
    #   plot = plot + 
    #     annotate("text", x = min(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "quadratic: NA", hjust=0, vjust = 2, color="red")
    # }
    # #print and plot tipping point/critical rate for positive change if possible
    # if(neg_rates_ss) { #if we can estimate a tipping point
    #   plot = plot + 
    #     annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = paste("spline:", round(tp_pt_rate,3)), hjust=1, vjust = 4, color="blue") +
    #     geom_point(aes(x=tp_pt_rate, y=tp_pt_growth), color="blue", size = 3) 
    # } else {
    #   plot = plot + 
    #     annotate("text", x = max(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "spline: NA", hjust=1, vjust = 4, color="blue")
    # }
    # #print and plot tipping point/critical rate for negative change if possible
    # if(neg_rates_ss_min) { #if we can estimate a tipping point
    #   plot = plot + 
    #     annotate("text", x = min(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = paste("spline:", round(tp_pt_rate_min,3)), hjust=0,vjust = 4, color="blue") +
    #     geom_point(aes(x=tp_pt_rate_min, y=tp_pt_growth_min), color="blue", size = 3) 
    # } else {
    #   plot = plot + 
    #     annotate("text", x = min(c(evo_scale * dY_quad, evo_scale * dY_ss)), y = max(c(input$rmax - (quad_max - pred_quad)/input$gT, input$rmax - (ss_max - pred_ss)/input$gT)), label = "spline: NA", hjust=0, vjust = 4, color="blue")
    # }
    # plot
    
    x <- with(data, seq(min(phenotype), max(phenotype), (max(phenotype)-min(phenotype))/params$nx)) #divide up the x-axis
    
    #constant to get rate of evo from selction gradient
    evo_scale = params$V_a / ((params$B - params$b) * params$gT) 
    
    #refit on x-centred data (to align with derivative)
    pred_quad <- log(quadratic.estimator(data,rowMeans(embed(x,2)))) #quad
    quad_max <- max(pred_quad) #max fitness
    pred_spline <- log(spline.estimator(data,rowMeans(embed(x,2)))) #spline
    spline_max <- max(pred_spline) #max fitness

    #perform bootstrap
    q.d.cis <- as.data.frame(quad.deriv.cis(data, x, B=params$bs, alpha=0.05))
    sp.d.cis <- as.data.frame(spline.deriv.cis(data, x, B=params$bs, alpha=0.05))    

    #get optimal trait values
    quad_opt <- x[which.max(quadratic.estimator(data,x))]
    spline_opt <- x[which.max(spline.estimator(data,x))]
    
    #maxima in derivative over positive lags
    pos_lags <- which(spline_opt - sp.d.cis$x > 0) #index for positive lags
    pos_dY_ss <- sp.d.cis$main.curve[pos_lags] #derivatives over positive lags
    pos_maxs <- which(diff(sign(diff(pos_dY_ss)))==-2) + 1 #local maxima in derivative over positive lags
    pos_nmaxs <- length(pos_maxs) #number of maxima
    
    #minima in derivative over negative lags
    neg_lags <- which(spline_opt - sp.d.cis$x < 0) #index for negative lags
    neg_dY_ss <- sp.d.cis$main.curve[neg_lags] #derivatives over negative lags
    neg_mins <- which(diff(sign(diff(neg_dY_ss)))==2) + 1 #local minima in derivative over negative lags
    neg_nmins <- length(neg_mins) #number of minima
        
    #find critical rates, if possible, from quadratic
    #first get indices for positive growth rates
    pos_grow_quad <- which(params$rmax - (quad_max - pred_quad) / params$gT > 0) #list of rows that have positive growth
    #look at positive rates of change
    crit_rate <- max(evo_scale * q.d.cis$main.curve[pos_grow_quad]) #potential critical rate (max rate of evo with positive growth)
    crit_index <- which.max(evo_scale * q.d.cis$main.curve[pos_grow_quad]) #get the index of the potential critical rate
    crit_growth <- params$rmax - (quad_max - pred_quad[pos_grow_quad][crit_index]) / params$gT #get the growth rate at the potential critical rate
    neg_rates_quad <- any(evo_scale * q.d.cis$main.curve > 0 & params$rmax - (quad_max - pred_quad) / params$gT < 0) #are there positive rates of change with negative growth rates (ie can we esimtate a critical rate or does the data not extend far enough)
    #look at negative rates of change
    crit_rate_min <- min(evo_scale * q.d.cis$main.curve[pos_grow_quad]) #potential critical rate (min rate of evo with positive growth)
    crit_index_min <- which.min(evo_scale * q.d.cis$main.curve[pos_grow_quad]) #get the index of the potential critical rate
    crit_growth_min <- params$rmax - (quad_max - pred_quad[pos_grow_quad][crit_index_min]) / params$gT #get the growth rate at the potential critical rate
    neg_rates_quad_min <- any(evo_scale * q.d.cis$main.curve < 0 & params$rmax - (quad_max - pred_quad) / params$gT < 0) #are there negative rates of evo with negative growth rates (ie can we esimtate a critical rate or does the data not extend far enough)
    
    #find critical rates/tipping points from spline
    #pos growth rates
    pos_grow <- which(params$rmax - (spline_max - pred_spline) / params$gT > 0) #list of rows that have positive growth
    #positive rates of change
    tp_pt_rate <- max(evo_scale * sp.d.cis$main.curve[pos_grow]) #potential tipping point
    tp_index <- which.max(evo_scale * sp.d.cis$main.curve[pos_grow]) #get the index
    tp_pt_growth <- params$rmax - (spline_max - pred_spline[pos_grow][tp_index]) / params$gT #growth rate at tipping point
    neg_rates_ss <- any(evo_scale * sp.d.cis$main.curve > 0 & params$rmax - (spline_max - pred_spline) / params$gT < 0) #are there positive rates of change with negative growth rate (ie can we esimtate a tipping point)
    #negative rates of change
    tp_pt_rate_min <- min(evo_scale * sp.d.cis$main.curve[pos_grow]) #potential tipping point
    tp_index_min <- which.min(evo_scale * sp.d.cis$main.curve[pos_grow]) #get the index
    tp_pt_growth_min <- params$rmax - (spline_max - pred_spline[pos_grow][tp_index_min]) / params$gT #growth rate at tipping point
    neg_rates_ss_min <- any(evo_scale * sp.d.cis$main.curve < 0 & params$rmax - (spline_max - pred_spline) / params$gT < 0, na.rm=T) #are there negative rates of evo with negative growth rates (ie can we esimtate a tipping point)
    
    #separate the positive and negative rates of change so that we have two functions for each fit, allowing geom_ribbon to work (requires a function, like geom_line)
    pos.q.rates <- which(evo_scale * q.d.cis$main.curve >= 0)
    neg.q.rates <- which(evo_scale * q.d.cis$main.curve < 0)
    pos.sp.rates <- which(evo_scale * sp.d.cis$main.curve >= 0)
    neg.sp.rates <- which(evo_scale * sp.d.cis$main.curve < 0)
    
    #rate of popn growth as a function of the rate of environmental change
    xmin = params$min_growth
    xmax = params$max_growth
    ymin = params$min_rate
    ymax = params$max_rate
    plot =
      ggplot() +
      geom_path(data = q.d.cis, aes(x = params$rmax - (quad_max - pred_quad)/params$gT, y = evo_scale * main.curve), colour="red", size=1) +
      geom_ribbon(aes(x = params$rmax - (quad_max - pred_quad[pos.q.rates])/params$gT, ymin = evo_scale * q.d.cis$lower.ci[pos.q.rates], ymax = evo_scale * q.d.cis$upper.ci[pos.q.rates]), fill="red", alpha="0.5") +
      geom_ribbon(aes(x = params$rmax - (quad_max - pred_quad[neg.q.rates])/params$gT, ymin = evo_scale * q.d.cis$lower.ci[neg.q.rates], ymax = evo_scale * q.d.cis$upper.ci[neg.q.rates]), fill="red", alpha="0.5") +
      # geom_path(data = sp.d.cis, aes(x = input$rmax - (spline_max - pred_spline)/input$gT, y = evo_scale * main.curve), colour="blue") +
      geom_ribbon(aes(x = params$rmax - (spline_max - pred_spline[pos.sp.rates])/params$gT, ymin = evo_scale * sp.d.cis$lower.ci[pos.sp.rates], ymax = evo_scale * sp.d.cis$upper.ci[pos.sp.rates]), fill="blue", alpha="0.5") +
      geom_ribbon(aes(x = params$rmax - (spline_max - pred_spline[neg.sp.rates])/params$gT, ymin = evo_scale * sp.d.cis$lower.ci[neg.sp.rates], ymax = evo_scale * sp.d.cis$upper.ci[neg.sp.rates]), fill="blue", alpha="0.5") +
      labs(x = 'long-term mean Malthusian growth rate', y = 'rate of environmental change') +
      theme_bw() + 
      theme(panel.border= element_blank()) +
      geom_hline(yintercept = 0) + 
      geom_vline(xintercept = 0) +
      xlim(xmin, xmax) + 
      ylim(ymin, ymax) +
      annotate("text", x = xmax, y = ymax, label = "~underline('max sustainable rates')", parse = TRUE, hjust=1, vjust = 0, color="black") +
      annotate("text", x = xmax, y = ymin, label = "~underline('min sustainable rates')", parse = TRUE, hjust=0, vjust = 0, color="black") +
      coord_flip()
    #plot spline over positive lags
    if(pos_nmaxs == 0) {
      #if no maxs then just plot solid curve over all lags
      plot = plot + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[pos_lags])/params$gT, y = evo_scale * sp.d.cis$main.curve[pos_lags]), colour="blue", size=1)
    } else {
      #if some maxs first plot dashed curve over all lags and solid curve from lag of zero to first max
      plot = plot + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[pos_lags])/params$gT, y = evo_scale * sp.d.cis$main.curve[pos_lags]), colour="blue", lty=2) + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[pos_lags][length(pos_lags):pos_maxs[pos_nmaxs]])/params$gT, y = evo_scale * sp.d.cis$main.curve[pos_lags][length(pos_lags):pos_maxs[pos_nmaxs]]), colour="blue", size=1) 
      keep_going <- TRUE
      current_index <- pos_maxs[pos_nmaxs]
      while(keep_going) {
        if(any(sp.d.cis$main.curve[pos_lags][current_index] < sp.d.cis$main.curve[pos_lags] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags])) {
          #if there is a greater lag that gives a greater derivative, set this to be the next stable lag
          next_index <- max(which(sp.d.cis$main.curve[pos_lags][current_index] < sp.d.cis$main.curve[pos_lags] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags]))
          if(any(sp.d.cis$main.curve[pos_lags][current_index] < sp.d.cis$main.curve[pos_lags][pos_maxs] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags][pos_maxs])) {
            #if there is another max at a greater lag, plot the solid line only up to there
            last_index <- pos_maxs[max(which(sp.d.cis$main.curve[pos_lags][current_index] < sp.d.cis$main.curve[pos_lags][pos_maxs] & spline_opt - sp.d.cis$x[pos_lags][current_index] < spline_opt - sp.d.cis$x[pos_lags][pos_maxs]))]
            current_index <- last_index
          } else {
            #if there are no more maxs then plot all the way to the end of the data and stop the loop
            last_index <- 1
            keep_going <- FALSE
          }
          #make temporary data frame with rate of enviro change and steady-state lag over range of interest (this lets ggplot know we want to plot a new layer each time)
          gg.data <- data.frame(x = params$rmax - (spline_max - pred_spline[pos_lags][next_index:last_index])/params$gT, y = evo_scale * sp.d.cis$main.curve[pos_lags][next_index:last_index])
          plot = plot +
            geom_path(data = gg.data, aes(x = x, y = y), colour="blue", size=1)
        } else {
          #if there are no more maxs then stop the loop
          keep_going <- FALSE
        }
      }  
    } 
    #plot spline over negative lags
    if(neg_nmins == 0) {
      #if no mins then just plot solid curve over all lags
      plot = plot + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[neg_lags])/params$gT, y = evo_scale * sp.d.cis$main.curve[neg_lags]), colour="blue", size=1)
    } else {
      #if some mins first plot dashed curve over all lags and solid curve to first min
      plot = plot + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[neg_lags])/params$gT, y = evo_scale * sp.d.cis$main.curve[neg_lags]), colour="blue", lty=2) + 
        geom_path(aes(x = params$rmax - (spline_max - pred_spline[neg_lags][1:neg_mins[1]])/params$gT, y = evo_scale * sp.d.cis$main.curve[neg_lags][1:neg_mins[1]]), colour="blue", size=1) 
      keep_going <- TRUE
      current_index <- neg_mins[1]
      while(keep_going) {
        if(any(sp.d.cis$main.curve[neg_lags][current_index] > sp.d.cis$main.curve[neg_lags] & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags])) {
          next_index <- min(which(sp.d.cis$main.curve[neg_lags][current_index] > sp.d.cis$main.curve[neg_lags] & spline_opt - sp.d.cis$x[neg_lags][current_index] > spline_opt - sp.d.cis$x[neg_lags]))
          if(any(sp.d.cis$main.curve[neg_lags][next_index] > sp.d.cis$main.curve[neg_lags][neg_mins] & spline_opt - sp.d.cis$x[neg_lags][next_index] > spline_opt - sp.d.cis$x[neg_lags][neg_mins])) {
            last_index <- neg_mins[min(which(sp.d.cis$main.curve[neg_lags][next_index] > sp.d.cis$main.curve[neg_lags][neg_mins] & spline_opt - sp.d.cis$x[neg_lags][next_index] > spline_opt - sp.d.cis$x[neg_lags][neg_mins]))]
            current_index <- last_index
          } else {
            last_index <- length(neg_lags)
            keep_going <- FALSE
          }
          gg.data <- data.frame(x = params$rmax - (spline_max - pred_spline[neg_lags][next_index:last_index])/params$gT, y = evo_scale * sp.d.cis$main.curve[neg_lags][next_index:last_index])
          plot = plot +
            geom_path(data = gg.data, aes(x = x, y = y), colour="blue", size=1)
        } else {
          keep_going <- FALSE
        }
      }  
    }
    #print and plot critical rate for positive change if possible
    if(neg_rates_quad) { #if we can estimate a critical rate
      plot = plot + 
        annotate("text", x = xmax, y = ymax, label = paste("quadratic:", round(crit_rate,3)), hjust=1, vjust = 2, color="red") +
        geom_point(aes(x = 0, y = crit_rate), color="red", size = 3) 
    } else {
      plot = plot + 
        annotate("text", x = xmax, y = ymax, label = "quadratic: NA", hjust=1, vjust = 2, color="red")
    }
    #print and plot critical rate for negative change if possible
    if(neg_rates_quad_min) { #if we can estimate a critical rate
      plot = plot + 
        annotate("text", x = xmax, y = ymin, label = paste("quadratic:", round(crit_rate_min,3)), hjust=0, vjust = 2, color="red") +
        geom_point(aes(x = 0, y = crit_rate_min), color="red", size = 3) 
    } else {
      plot = plot + 
        annotate("text", x = xmax, y = ymin, label = "quadratic: NA", hjust=0, vjust = 2, color="red")
    }
    #print and plot tipping point/critical rate for positive change if possible
    if(neg_rates_ss) { #if we can estimate a tipping point
      plot = plot + 
        annotate("text", x = xmax, y = ymax, label = paste("spline:", round(tp_pt_rate,3)), hjust=1, vjust = 4, color="blue") +
        geom_point(aes(x = tp_pt_growth, y = tp_pt_rate), color="blue", size = 3) 
    } else {
      plot = plot + 
        annotate("text", x = xmax, y = ymax, label = "spline: NA", hjust=1, vjust = 4, color="blue")
    }
    #print and plot tipping point/critical rate for negative change if possible
    if(neg_rates_ss_min) { #if we can estimate a tipping point
      plot = plot + 
        annotate("text", x = xmax, y = ymin, label = paste("spline:", round(tp_pt_rate_min,3)), hjust=0,vjust = 4, color="blue") +
        geom_point(aes(x = tp_pt_growth_min, y = tp_pt_rate_min), color="blue", size = 3) 
    } else {
      plot = plot + 
        annotate("text", x = xmax, y = ymin, label = "spline: NA", hjust = 0, vjust = 4, color="blue")
    }
    plot
    
  })  
}

# Run the app ----
shinyApp(ui = ui, server = server)
