#---Import Dependencies
library(dplyr)
library(ggplot2)
library(ggforce)
library(geospark)
library(ggmap)
library(e1071)
library(MLmetrics)
library(gridExtra)
library(randomForest)
library(caret)
library(tidyr)

#----Preprocessing----
setwd()
data <- read.csv("data.csv")

#assign region order
data$Region <- factor(data$Region, levels=c("NWIP","DML","MSP","MOP", "SIDP","IS"))
regions <- levels(data$Region)

#filter out null predictor and nitrates + convert warning to factor
data <- drop_na(data)
data$Warning_level <- factor(data$Warning_level, levels = c("Safe","Warning","Alert"))
data$Future_Warning <- factor(data$Future_Warning, levels = c("Safe","Warning","Alert"))

#create train/test splitting function
define_train_test <- function(){
  train <- list()
  test <- list()
  
  #assign train and test set based on year
  for (i in 1:6){
    train[[i]] <- data[data$Year != 2025 & data$Region == regions[i],]
    
    test[[i]] <- data[data$Year == 2025 & data$Region == regions[i],]
  }
  return (list(train,test))
}

#create train and test variables
train_test <- define_train_test()
train <- train_test[[1]]
test <- train_test[[2]]

#---Use Cross Validation to determine SVM hyperparameters
# Initialize CV object
cv_svm_models_rbf <- list()

#set seed
set.seed(12345)

# Define training control with 5-fold cross-validation
train_control <- trainControl(
  method = "cv",       # Cross-validation
  number = 5,         # Number of folds
  verboseIter = TRUE   # Show progress
)

# Define parameter grid for tuning
svm_grid_rbf <- expand.grid(
  C = c(0.1, 1, 10, 100),
  sigma = c(0.01, 0.1, 1)  # caret uses 'sigma' instead of 'gamma'
)

#find the best parameters via CV
for (i in 1:6){
  # Train SVM with radial basis function (RBF) kernel
  cv_svm_models_rbf[[i]] <- train(
    Change ~ Monthly_Rainfall+ Avg_temp, 
    data = train[[i]],
    method = "svmRadial",
    trControl = train_control,
    tuneGrid = svm_grid_rbf,
    preProcess = c("center", "scale"), 
    metric = "RMSE"
  )
  
}

# Define the best parameters for the first iteration
svm_params_rbf <- cv_svm_models_rbf[[1]]$bestTune

# Add to this list for the remaining items
for (i in 2:length(cv_svm_models_rbf)){
  svm_params_rbf <- rbind(svm_params_rbf, cv_svm_models_rbf[[i]]$bestTune)
}

#---Training the Models---
#initialize lists to store model objects
lm_models <- list()
svm_models_rbf <- list()
lm_predicted_change <- list()
svm_predicted_change <- list()

#train polynomial model via bootstrap
for (i in 1:6){
  model <- list()
  pred <- list()
  set.seed(35791)
  
  for (j in 1:1000){
    indices <- sample(nrow(train[[i]]), nrow(train[[i]]), replace=TRUE)
    model[[j]] <- lm(Change ~ 0+poly(Monthly_Rainfall, degree=2) + poly(Avg_temp, degree=2), data=train[[i]][indices,])
    
    pred[[j]] <- predict(model[[j]],test[[i]])
  }
  
  transpose <- t(sapply(pred,function(x){
    x
  }))
  
  lm_predicted_change[[i]] <- colMeans(transpose)
  
}

#train svm model - rbf via bootstrap
for (i in 1:6){
  model <- list()
  pred <- list()
  set.seed(35791)
  
  for (j in 1:1000){
    #sample indices with replacement
    indices <- sample(nrow(train[[i]]), nrow(train[[i]]), replace=TRUE)
    
    #run model for all j 
    model[[j]] <- svm(Change ~ 0+Monthly_Rainfall+Avg_temp, data = train[[i]][indices,], kernel = "radial",
                               cost = svm_params_rbf$C[i], gamma = svm_params_rbf$sigma[i], scale = TRUE,
                               type = 'eps-regression')
    
    #predict
    pred[[j]] <- predict(model[[j]],test[[i]])
  }
  
  #average all predictions for each model
  transpose <- t(sapply(pred,function(x){
    x
  }))
  
  svm_predicted_change[[i]] <- colMeans(transpose)

}

#---Evaluate Model Performance
get_evaluators <- function(model, predicted_change){
  pred_warning <- list() #predicted warning level (Safe, Warning, Alert)
  predicted_outputs <- list() #expected nitrate concentration (mg/L) -- change + prev month
  resids <- list() #residuals of expected - actual change
  f1s <- list()
  recalls <- list()
  accuracy <- c()
  
  #iterate over each region
  for (i in 1:6){
    #Note: Remove nulls, which have been filtered to only be missing change variable, meaning
    #a prediction would not be able to be compared to the row, thus not used for testing. 
    #compute actual nitrate output
    predicted_outputs[[i]] <- test[[i]]$Avg_Nitrate + predicted_change[[i]] 
    
    #convert to warning class --> some low instances ended up less than zero, thus defaulted to safe class with min as the lower endpoint
    pred_warning[[i]] <- cut(predicted_outputs[[i]],breaks = c(-Inf,7.5,10,Inf), labels = c("Safe","Warning","Alert"),right = FALSE)
    
    #compute residuals
    resids[[i]] <- test[[i]]$Change - predicted_change[[i]]
    
    #f1s
    f1_low <- F1_Score(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                       positive = "Safe")
    f1_med <- F1_Score(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                       positive = "Warning")
    f1_high <- F1_Score(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                        positive = "Alert")
    
    #recalls
    recall_low <- Recall(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                         positive = "Safe")
    recall_med <- Recall(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                         positive = "Warning")
    recall_high <- Recall(y_true = test[[i]]$Future_Warning, y_pred = pred_warning[[i]],
                          positive = "Alert")
    
    #accuracy
    accuracy[i] <- Accuracy(y_true = drop_na(test[[i]])$Future_Warning, y_pred = pred_warning[[i]])
    
    f1s[[i]] <- c(f1_low, f1_med, f1_high)
    recalls[[i]] <- c(recall_low, recall_med, recall_high)
    
  }
  return (list(f1s, recalls, accuracy, predicted_outputs, pred_warning, resids))
}

#initialize evaluators for linear model
linear_evals <- get_evaluators(lm_models, lm_predicted_change)
linear_f1s <- linear_evals[[1]]
linear_recalls <- linear_evals[[2]]
linear_accuracy <- linear_evals[[3]]
linear_predicted <- linear_evals[[4]]
linear_pred_warning <- linear_evals[[5]]
linear_resids <- linear_evals[[6]]

#initialize evaluators for SVM model
svm_mod_rbf <- get_evaluators(svm_models_rbf, svm_predicted_change)
svm_f1s_rbf <- svm_mod_rbf[[1]]
svm_recalls_rbf <- svm_mod_rbf[[2]]
svm_accuracy_rbf <- svm_mod_rbf[[3]]
svm_predicted_rbf <- svm_mod_rbf[[4]]
svm_pred_warning_rbf <- svm_mod_rbf[[5]]
svm_residuals_rbf <- svm_mod_rbf[[6]]

#fill in f1 nulls to 0 
for (i in 1:6){
  for (j in 1:3){
    if (is.na(linear_f1s[[i]][j])){
      linear_f1s[[i]][j] <- 0
    }
    if (is.na(svm_f1s_rbf[[i]][j])){
      svm_f1s_rbf[[i]][j] <- 0
    }
    
  }
}

#---Plot Results
plot_models <- function(test, predicted_outputs_list, resids, model_type){
  region_dfs <- list()
  conf_matrix <- list()
  plots <- list()
  res_plots <- list()
  
  for (i in 1:6){
    #df of actual against predicted and warning levels
    region_dfs[[i]] <- data.frame(Actual = drop_na(test[[i]])$Future_Nitrate,
                                  Predicted = predicted_outputs_list[[i]], 
                                  Actual_warning = drop_na(test[[i]])$Future_Warning,
                                  Predicted_warning = cut(as.numeric(predicted_outputs_list[[i]]), 
                                                          breaks = c(min(predicted_outputs_list[[i]])-1,7.5,10,30),
                                                          labels = c("Safe","Warning","Alert"),right=FALSE),
                                  Residuals = resids[[i]])
    
    #confusion matrices
    conf_matrix[[i]] <- table(Actual = region_dfs[[i]]$Actual_warning,
                              Predicted = region_dfs[[i]]$Predicted_warning)
    
    #plots of pred v. actual
    plots[[i]] <- ggplot(region_dfs[[i]])+
      geom_point(aes(x=Predicted,y=Actual), color = "darkred")+
      labs(title = paste0("Region: ",regions[i], " - ", model_type))+
      geom_abline(slope=1, intercept = 0, linetype="dashed", color="darkblue")+
      theme(plot.title = element_text(hjust=.5, face="bold"))
    
    #plots of residuals v. actual
    res_plots[[i]] <- ggplot(region_dfs[[i]])+
      geom_point(aes(x=Actual, y = Residuals),color="darkgreen")+
      labs(x="Actual Value",y="Residual",title = paste0("Region: ",regions[i], " Residuals"))+
      geom_abline(slope=0, intercept = 0, linetype="dashed", color="darkblue")+
      theme(plot.title = element_text(hjust=.5, face="bold"))
    
  }
  return (list(region_dfs,conf_matrix, plots, res_plots))
}

#initialize plot elements for lm
linear_plot_func <- plot_models(test, linear_predicted, linear_resids, "Polynomial")
linear_results <- linear_plot_func[[1]]
linear_cm <- linear_plot_func[[2]]
linear_plot <- linear_plot_func[[3]]
linear_resid_plot <- linear_plot_func[[4]]

#initialize plot elements for svm
svm_plot_func <- plot_models(test, svm_predicted_rbf, svm_residuals_rbf, "SVR")
svm_results <- svm_plot_func[[1]]
svm_cm <- svm_plot_func[[2]]
svm_plots <- svm_plot_func[[3]]
svm_resid_plots <- svm_plot_func[[4]]

#---Compare models via Wilcoxon Signed Ranks
calculate_wilcox <- function(){
  #--F1s
  f1_low_lm <- unlist(linear_f1s)[seq(from=1, to = length(unlist(linear_f1s)),by=3)]
  f1_med_lm <- unlist(linear_f1s)[seq(from=2, to = length(unlist(linear_f1s)),by=3)]
  f1_high_lm <- unlist(linear_f1s)[seq(from=3, to = length(unlist(linear_f1s)),by=3)]
  
  f1_low_svm <- unlist(svm_f1s_rbf)[seq(from=1, to = length(unlist(svm_f1s_rbf)),by=3)]
  f1_med_svm <- unlist(svm_f1s_rbf)[seq(from=2, to = length(unlist(svm_f1s_rbf)),by=3)]
  f1_high_svm <- unlist(svm_f1s_rbf)[seq(from=3, to = length(unlist(svm_f1s_rbf)),by=3)]
  
  wilcox_f1_low <- wilcox.test(f1_low_lm, f1_low_svm, paired = TRUE, alternative= "two.sided")
  wilcox_f1_med <- wilcox.test(f1_med_lm, f1_med_svm, paired = TRUE, alternative= "two.sided")
  wilcox_f1_high <- wilcox.test(f1_high_lm, f1_high_svm, paired = TRUE, alternative= "two.sided")
  
  #--Accuracy
  wilcox_acc <- wilcox.test(linear_accuracy, svm_accuracy_rbf, paired = TRUE, alternative = "two.sided")
  
  #--Recalls
  recall_low_lm <- unlist(linear_recalls)[seq(from=1, to = length(unlist(linear_recalls)),by=3)]
  recall_med_lm <- unlist(linear_recalls)[seq(from=2, to = length(unlist(linear_recalls)),by=3)]
  recall_high_lm <- unlist(linear_recalls)[seq(from=3, to = length(unlist(linear_recalls)),by=3)]
  
  recall_low_svm <- unlist(svm_recalls_rbf)[seq(from=1, to = length(unlist(svm_recalls_rbf)),by=3)]
  recall_med_svm <- unlist(svm_recalls_rbf)[seq(from=2, to = length(unlist(svm_recalls_rbf)),by=3)]
  recall_high_svm <- unlist(svm_recalls_rbf)[seq(from=3, to = length(unlist(svm_recalls_rbf)),by=3)]
  
  wilcox_recall_low <- wilcox.test(recall_low_lm, recall_low_svm, paired = TRUE, alternative= "two.sided")
  wilcox_recall_med <- wilcox.test(recall_med_lm, recall_med_svm, paired = TRUE, alternative= "two.sided")
  wilcox_recall_high <- wilcox.test(recall_high_lm, recall_high_svm, paired = TRUE, alternative= "two.sided")
  
  return (list(wilcox_f1_low, wilcox_f1_med, wilcox_f1_high, wilcox_acc, wilcox_recall_low, wilcox_recall_med, wilcox_recall_high))
}

wilcox <- calculate_wilcox()

#get r_squared for goodness of fit
lm_r <- lm(Predicted~Actual,data=linear_results[[2]])
svm_r <- lm(Predicted~Actual, data=svm_results[[2]])

#add R^2 to plots for DML region
linear_plot[[2]] <- linear_plot[[2]] + annotate("text",x=3, y=25, label = paste0("R-squared: ",round(summary(lm_r)$r.squared,3)))
svm_plots[[2]] <- svm_plots[[2]] + annotate("text",x=4, y=25, label = paste0("R-squared: ",round(summary(svm_r)$r.squared,3)))

grid.arrange(linear_plot[[2]],svm_plots[[2]],ncol=2)
