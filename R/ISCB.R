
### Author: Saed Sayad M.D. Ph.D.; saed@bioada.com ###

# Load libraries
library(data.table)
library(formattable)
library(plotrix)
library(limma)
library(dplyr)
library(Rtsne)
library(MASS)
library(xgboost)
library(rafalib)
library(factoextra)
library(caTools)
library(pROC)
library(caret)
library(gains)
library(lift)

#Remove ALL
rm(list=ls())

####################
# Data Preparation #
####################


# Read Expressions file
df <- read.csv("C:\\Temp\\ISCB\\GSE74763_rawlog_expr.csv")
df2 <- df[,-1]
rownames(df2) <- df[,1]
expr <- transpose(df2)
rownames(expr) <- colnames(df2)
colnames(expr) <- rownames(df2)
dim(expr)

# Read Samples file
targets <- read.csv("C:\\Temp\\ISCB\\GSE74763_rawlog_targets.csv")
colnames(targets)
dim(targets)

# Merge Expressions with Samples
data <- cbind(expr, targets)
colnames(data)
dim(data)


#######################
# Univariate Anlaysis #
#######################

# Density plot for all expressions
plotDensities(expr, main="Expressions", legend=F)

# Categorical Variables -  Count and Count%
univar.count <- table(data$target)/nrow(data)
univar.countpct <- formattable::percent(univar.count)
univar.count
univar.countpct

# Categorical Variables - Pie and Bar charts
pie(table(data$target), main = "target")
lbl <- unique(data$target)
pie3D(table(data$target),labels=lbl, explode = 0.1, main = "target")
barplot(table(data$target), main="target", xlab="", ylab="Count", col="blue")

# Numerical Variables - Summary (Descriptive Statistics)
summary(data$age)

# Numerical Variables - Plots
boxplot(data$age, col="cyan", main="age", ylab="intensity")
hist(data$age, col="cyan", main="age", ylab="intensity")
plotDensities(data$age, main="age", col="blue")


######################
# Bivariate Anlaysis #
######################

## Numerical & Numerical ##

# Linear Correlation
cor(data$P000001,data$P000002)

# Scatter plot
pairs(~P000001+P000002, data=data, main="GSE74763_rawlog", col="darkgreen")

## Categorical & Categorical ##

# Count (crosstab)
xtabs(~target+gender, data=data)

# CrossTab Plot
xt <- xtabs(~target+gender, data=data)
plot(xt, main="GSE74763_rawlog", col="darkgreen")

# Chi2 test
tbl <- table(data$target, data$gender)
chisq.test(tbl)

## Categorical & Numerical ##

# Box plot
boxplot(P000001~target, data=data, col="darkgreen", main="GSE74763_rawlog", ylab="P000001 - intensity")

# t-test
t.test(data$P000002~data$target)

# ANOVA
fit <- aov(P000001~target, data=data)
summary(fit)


#########################
# Multivariate Anlaysis #
#########################

# PCA (Principal Components Analysis)
d1 <- data
d1$geo_accession <- NULL
d1$gage <- NULL
d1$gender <- NULL
d1$target <- NULL

pc <- prcomp(d1)

groups <- as.fumeric(data$target)
plot(pc$x[, 1], pc$x[, 2], main = "PCA", col=groups, xlab = "PC1", ylab = "PC2")

# SVD
cx <- sweep(d1, 2, colMeans(d1), "-")
sv <- svd(cx)
names(sv)

plot(sv$d^2/sum(sv$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained", main="GSE74763_rawlog")

# t-SNE (t-distributed Stochastic Neighbour Embedding)
d1 <- data
d1$geo_accession <- NULL
d1$age <- NULL
d1$gender <- NULL
d1$target <- NULL
groups <- as.fumeric(data$target)

colors = rainbow(length(unique(groups)))
names(colors) = unique(groups)
tsne <- Rtsne(d1, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
plot(tsne$Y, col=colors[groups], main="GSE74763_rawlog")


#######################
# Predictive Modeling #
#######################

# Variable Selection using t-test to select the topN variables (probes) from 9480 probes.
tval <- lapply(d1, function(x) { t.test(x ~ groups)$statistic })

tdf <- as.data.frame(t(as.data.frame(tval)))
colnames(tdf)
tdf$t <- abs(tdf$t)
tdfo <- tdf[order(tdf$t,decreasing=TRUE),,drop=FALSE]
head(tdfo)

# Logistic Regression
d1$target <- as.fumeric(data$target)-1
model <- glm(target ~ P000833+P007414+P002449, data = d1, family = binomial, maxit = 100)
print(model)

# XGBoost (Extreme Gradient Boosting)
bst <- xgboost(data = as.matrix(d1), label = d1$target, max.depth = 5, eta = 1, nthread = 2, nrounds = 100
               , objective = "binary:logistic")

var_importance <- xgb.importance(model = bst)
print(var_importance)


######################
# Model Evaluation 1 #
######################

# Splitting Data into Training and Test Sets
set.seed(101)

sample = sample.split(data$target, SplitRatio = .8)
train = subset(data, sample == TRUE)
test = subset(data, sample == FALSE)

dim(train)
dim(test)

#
tr <- train
tr$geo_accession <- NULL
tr$age <- NULL
tr$gender <- NULL
groups <- as.fumeric(tr$target)
tr$target <- NULL


# Variable Selection using t-test to select the topN variables (probes) from 9480 probes.
tval <- lapply(tr, function(x) { t.test(x ~ groups)$statistic })

tdf <- as.data.frame(t(as.data.frame(tval)))
colnames(tdf)
tdf$t <- abs(tdf$t)
tdfo <- tdf[order(tdf$t,decreasing=TRUE),,drop=FALSE]
head(tdfo)


# Logistic Regression

#train
train$target <- as.fumeric(train$target)-1
model <- glm(target ~ P000833+P002449+P009369, data = train, family=binomial(logit), maxit = 100)
print(model)

#test
test$geo_accession <- NULL
test$age <- NULL
test$gender <- NULL
test$target <- as.fumeric(test$target)-1

pb <- predict(model, test, type="response")
pb <- as.data.frame(pb)

# Confusion Matrix
pc <- NULL
pc <- ifelse(pb$pb > 0.5,"1","0")
summary(pc)
xtab <- table(pc, test$target)
caret::confusionMatrix(xtab, positive = "1")

# ROC Chart
pb <- NULL
pb <- predict(model, test, type="response")
pb <- as.data.frame(pb)
labels <- test$target
scores <- pb$pb

plot(roc(labels, scores, direction="<"), col="blue", lwd=3, main="ROC Chart")

auc(roc(labels, scores, direction="<"))


# Gain and Lift Charts
pb <- NULL
pb <- predict(model, test, type="response")
pb <- as.data.frame(pb)
labels <- test$target
scores <- pb$pb

gains(labels, scores, groups=10)

plot(gains(labels, scores, groups=10))

plotLift(scores, labels, cumulative = TRUE, n.buckets = 10)


######################
# Model Evaluation 2 #
######################

# Splitting Data into Training and Test Sets
set.seed(101)

sample = sample.split(data$target, SplitRatio = .8)
train = subset(data, sample == TRUE)
test = subset(data, sample == FALSE)

dim(train)
dim(test)

#
tr <- train
tr$geo_accession <- NULL
tr$age <- NULL
tr$gender <- NULL
tr$target <- as.fumeric(tr$target)-1


# XGBoost (Extreme Gradient Boosting)
bst <- xgboost(data = as.matrix(tr), label = tr$target, max.depth = 5, eta = 1, nthread = 2, nrounds = 100
               , objective = "binary:logistic")

var_importance <- xgb.importance(model = bst)
print(var_importance)

# test
test$geo_accession <- NULL
test$age <- NULL
test$gender <- NULL
test$target <- as.fumeric(test$target)-1

pb <- predict(bst,as.matrix(test),missing='')

# Confusion Matrix
pc <- NULL
pc <- ifelse(pb > 0.5,"1","0")
summary(pc)
xtab <- table(pc, test$target)
caret::confusionMatrix(xtab, positive = "1")

# ROC Chart
pb <- NULL
pb <- predict(bst,as.matrix(test),missing='')
pb <- as.data.frame(pb)
labels <- test$target
scores <- pb$pb

plot(roc(labels, scores, direction="<"), col="blue", lwd=3, main="ROC Chart")

auc(roc(labels, scores, direction="<"))


# Gain and Lift Charts
pb <- NULL
pb <- predict(bst,as.matrix(test),missing='')
labels <- test$target
scores <- pb

gains(labels, scores, groups=10)

plot(gains(labels, scores, groups=10))

plotLift(scores, labels, cumulative = TRUE, n.buckets = 10)


