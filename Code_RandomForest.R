#construct the random forest model based on the GSE164174 
library(randomForest)

# loading dataset
data <- read.csv("inputfile_164174.csv",row.names = 1)

# training set: 70%; test set: 30%
split_data <- read.csv("164174_split.csv")
train_df = data[split_data$Train ,]
test_df = data[split_data$Test ,]
test_df <- na.omit(test_df)
train_df <- na.omit(train_df)
qpcr_df = read.csv("inputfile_qpcr.csv",row.names = 1) 

train_df$group = factor(train_df$group,
                        levels = c("0" , "1"))
test_df$group = factor(test_df$group ,
                       levels = c("0" , "1"))
qpcr_df$group = factor(qpcr_df$group ,
                       levels = c("0" , "1"))

set.seed(2345)
# select the best mtry parameters
tuneRF(train_df[,-1],train_df[,1],
       stepFactor = 0.5,
       plot = TRUE,
       ntreeTry = 1000,
       trace = TRUE,
       improve = 0.05)


set.seed(5678)
train_randomforest <- randomForest(
  group ~ .,  data = train_df,
  importance=TRUE ,
  proximity=TRUE)
train_randomforest

library(pROC)
fitted.prob<-predict(train_randomforest,
                     newdata=qpcr_df,type = "prob")
roc_prob<-roc(qpcr_df$group,
              fitted.prob[,2],ci=T,
              levels = c(0,1),
              direction="<") #ROC
auc(roc_prob)
roc_prob$ci
##confusion matrix
library(caret)
library(reportROC)

# prediction on the test
predicted_probs <- predict(train_randomforest,qpcr_df, type = "prob")  # 概率预测
predicted_classes <- predict(train_randomforest, qpcr_df, type = "response")  # 类别预测

#evaluate the performance
roc_results <- reportROC(
  gold = qpcr_df$group,  # true
  predictor = predicted_probs[, 2], 
  important = "se",
  plot = FALSE 
)

# print the results from reportROC
cat("reportROC Results:\n")
print(roc_results)

threshold <- 0.5
adjusted_predicted_classes <- ifelse(predicted_probs[, 2] >= threshold, "1", "0")

adjusted_predicted_classes <- factor(adjusted_predicted_classes, levels = c("0", "1"))
reference <- factor(qpcr_df$group,levels = c("0", "1"))

# confusion matrix
conf_matrix <- confusionMatrix(data = adjusted_predicted_classes, reference = reference,positive = "1")

cat("\nConfusion Matrix:\n")
print(conf_matrix)

# Sensitivity, Specificity
cat("\nMetrics:\n")
cat("Accuracy: ", conf_matrix$overall['Accuracy'], "\n")
cat("Sensitivity: ", conf_matrix$byClass['Sensitivity'], "\n")
cat("Specificity: ", conf_matrix$byClass['Specificity'], "\n")