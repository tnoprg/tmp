install.packages('rms')
install.packages('ggplot2')
install.packages('pROC')
install.packages('caret')
install.packages('brant')

library(MASS)        # 用于有序逻辑回归 (polr)
library(rms)         # 用于构建模型和列线图
library(ggplot2)     # 用于高级绘图（可选）
library(pROC)        # 用于ROC分析（可选）
library(caret)       # 用于数据分割（可选）
library(brant)       # 用于比例优势检验
library(Hmisc)  
library(dplyr)    # For data manipulation
library(boot) 
set.seed(123)

Rdata$丈夫学历 <- factor(Rdata$丈夫学历, 
                                 levels = c(1, 2, 3),
                                 labels = c("大专、本科", "研究生", "其他"))

summary(Rdata$丈夫学历)

Rdata$本人职业 <- factor(Rdata$本人职业, 
                     levels = c(1, 2, 3, 4, 5),
                   labels = c("体制工作","技术员","职员","个体户","其他"))

summary(Rdata$本人职业)
                   
Rdata$吸烟 <- factor(Rdata$吸烟, 
                   levels = c(1, 2),
                   labels  = c("是", "否") )
                  
summary(Rdata$吸烟)

Rdata$倾诉方式 <- factor(Rdata$倾诉方式, 
                     levels = c(1, 2, 3, 4),
                     labels = c("无", "被动", "1-2人", "主动") )

summary(Rdata$倾诉方式)   

Rdata$求助方式 <- factor(Rdata$求助方式, 
                     levels = c(1, 2, 3, 4),
                     labels = c("无", "很少", "有时", "经常")
                     )
summary(Rdata$求助方式)

Rdata$分娩次数 <- as.numeric(Rdata$分娩次数)
summary(Rdata$分娩次数)

Rdata$母乳 <- factor(Rdata$母乳,
                   levels = c(1, 2, 3),
                   labels = c("否","混合", "是")
                   )

summary(Rdata$母乳)

Rdata$抑郁等级 <- factor(Rdata$EPDS级别, 
                     levels = c(1, 2, 3, 4),
                     labels  = c("正常", "轻度", "中度", "重度"),
                     ordered = TRUE)
summary(Rdata$抑郁等级)


model <- polr(Rdata$抑郁等级 ~ SSRS + 
       PSQI  +  
       Rdata$分娩次数 + 
       TPOab +
       Rdata$孕酮 +
       Rdata$孕期EPDS +
       Rdata$丈夫学历 +
       Rdata$本人职业 +
       Rdata$吸烟 + 
       Rdata$倾诉方式 +
       Rdata$求助方式 +
       Rdata$母乳,
     data = Rdata,
     Hess = TRUE)

summary(model)

Ddist <- datadist(Rdata$SSRS, 
                  Rdata$PSQI,  
                    Rdata$分娩次数, 
                    Rdata$TPOab,
                    Rdata$孕酮,
                    Rdata$孕期EPDS,
                    Rdata$丈夫学历,
                    Rdata$本人职业,
                    Rdata$吸烟,
                    Rdata$倾诉方式,
                    Rdata$求助方式,
                    Rdata$母乳)  
options(datadist='Ddist')   
Rdata_clean <- na.omit(Rdata)

model1 <- lrm (抑郁等级 ~
                 SSRS + 
                 PSQI  +  
                 分娩次数 + 
                 TPOab +
                 孕酮 +
                 孕期EPDS +
                 丈夫学历 +
                 本人职业 +
                 吸烟 + 
                倾诉方式 +
                 求助方式 +
                母乳,data = Rdata_clean, x = FALSE, y = FALSE)



?(ctable <- coef(summary(model)))
p_values <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p_values))
print(ctable)


brant_test <- brant(model)
print(brant_test) 


predprob <- function(model, newdata) {
  preds <- predict(model, newdata = newdata, type = "probs")
  colMeans(preds)  # 返回平均预测概率
}

nom_pred <- function(model, newdata) {
  predict(model, newdata = newdata, type = "probs")
}

clean_data <- na.omit(Rdata)
summary(clean_data)
clean_data$抑郁等级 <- factor(clean_data$抑郁等级, ordered = TRUE)
dd <- datadist(clean_data)
options(datadist = "dd")
rms_model <- lrm(抑郁等级 ~ SSRS + PSQI +
                   丈夫学历+
                   本人职业+
                   吸烟+
                   倾诉方式+
                   求助方式+
                   分娩次数+
                   母乳+
                   TPOab+
                   孕酮+
                  孕期EPDS,
                 data = clean_data, x = TRUE, y = TRUE)
print(rms_model)

graph1 <- nomogram (rms_model, fun = plogis, lp = F,
                    funlabel = "产后抑郁风险")
plot(graph1)

validate_boot <- validate(rms_model, method = "boot", B = 1000)
print(validate_boot)

corrected_cindex <- validate_boot["Dxy", "index.corrected"]/2 + 0.5
cat("Optimism-corrected C-index:", corrected_cindex, "\n")

calibration <- calibrate(rms_model, method = "boot", B = 1000)
plot(calibration)

pred_probs <- predict(rms_model, type = "fitted")

depression_levels <- sort(unique(clean_data$抑郁等级))

binary_outcome <- as.numeric(clean_data$抑郁等级 > min(depression_levels))

roc_obj <- roc(binary_outcome, pred_probs)
auc_value <- auc(roc_obj)
ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 1000)
cat("AUC:", round(auc_value, 3), 
    "(95% CI:", round(ci[1], 3), "-", round(ci[3], 3), ")\n")

roc_list <- list()
auc_results <- data.frame(
  Threshold = character(),
  AUC = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  stringsAsFactors = FALSE
)
thresholds <- depression_levels[-length(depression_levels)]
summary(thresholds)

  for(i in seq_along(thresholds)) {
    # Convert to binary at this threshold
    binary_outcome <- as.numeric(clean_data$抑郁等级 >= thresholds[i])
    
    # Calculate ROC curve
    roc_obj <- roc(binary_outcome, pred_probs, quiet = TRUE)
    
    # Store ROC object
    roc_list[[i]] <- roc_obj
    
    # Calculate 95% CI for AUC using bootstrap
    set.seed(123 + i)  # For reproducibility
    ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)
    
    # Store results
    auc_results <- rbind(auc_results, data.frame(
      Threshold = paste("≥", thresholds[i]),
      AUC = as.numeric(auc(roc_obj)),
      Lower_CI = ci[1],
      Upper_CI = ci[3],
      stringsAsFactors = FALSE
    ))
  }  

print(auc_results)

simple_roc <- roc(as.numeric(clean_data$抑郁等级), pred_scores)
plot(simple_roc, main = "Simplified ROC Curve (Ordinal as Continuous)")
auc_simple <- auc(simple_roc)
ci_simple <- ci.auc(simple_roc)
cat("Simplified AUC:", round(auc_simple, 3), 
    "\n95% CI: [", round(ci_simple[2], 3), ",", round(ci_simple[3], 3), "]\n")


boot_auc <- function(data, indices) {
  # 从原始数据中抽样
  sampled_data <- data[indices, ]
  
  # 在抽样数据上重新拟合模型
  boot_model <- lrm(抑郁等级 ~ SSRS + PSQI + 分娩次数 + TPOab + 孕酮 + 孕期EPDS + 
                      丈夫学历 + 本人职业 + 吸烟 + 倾诉方式 + 求助方式 + 母乳, 
                    data = sampled_data, x = TRUE, y = TRUE)
  
  # 预测概率（在原始数据上预测，避免乐观偏差）
  pred_prob <- predict(boot_model, newdata = Rdata_clean, type = "fitted")
  
  # 计算ROC和AUC
  roc_obj <- roc(Rdata_clean$抑郁等级, pred_prob)
  return(as.numeric(roc_obj$auc))
}

set.seed(123)
boot_results <- boot(data = Rdata_clean, statistic = boot_auc, R = 1000)


  # 计算ROC和AUC
  roc_obj <- roc(Rdata_clean$抑郁等级, pred_prob)
  return(as.numeric(roc_obj$auc))
}

set.seed(123)
boot_results <- boot(data = sample_data, statistic = boot_auc, R = 1000)

mean_auc <- mean(boot_results$t)
auc_ci <- boot.ci(boot_results, type = "bca")
cat("Bootstrap AUC:", mean_auc, "\n")
cat("AUC 95% CI:", auc_ci$bca[4:5], "\n")

boot_auc <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # 重新拟合模型
  boot_model <- lrm(抑郁等级 ~ SSRS + PSQI + 分娩次数 + TPOab + 孕酮 + 孕期EPDS + 
                      丈夫学历 + 本人职业 + 吸烟 + 倾诉方式 + 求助方式 + 母乳, 
                    data = sampled_data, x = TRUE, y = TRUE)
  
  # 关键修正：在同一个 sampled_data 上计算预测概率和真实标签
  pred_probs <- predict(boot_model, newdata = sampled_data, type = "fitted")  # 预测概率
  binary_outcome <- as.numeric(sampled_data$抑郁等级) - 1  # 真实标签（转换为0/1）
  
  # 检查长度（此时应该完全一致）
  if (length(binary_outcome) != length(pred_probs)) {
    stop(paste("致命错误：长度不一致！binary_outcome:", length(binary_outcome), 
               "pred_probs:", length(pred_probs)))
  }
  
  # 计算ROC和AUC
  roc_obj <- roc(binary_outcome, pred_probs)
  return(as.numeric(roc_obj$auc))
}

dim(Rdata_clean)
set.seed(123)
sample_indices <- sample(1:nrow(Rdata_clean), replace = TRUE)
sampled_data <- Rdata_clean[sample_indices, ]
dim(sampled_data)
