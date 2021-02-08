########################################
# A general feature selection protocol #
# ==================================== #
# written by:                          #
#     Thi Ngoc Ha Nguyen               #
#     Haoliang Xue                     #
########################################

set.seed(91400)

library(glmnet)
library(randomForest)
library(parallel)

extractStableFeatures <- function(x, y, thres, num.runs, glmnet.alpha) {
    n <- nrow(x) 
    p <- ncol(x)
    features <- colnames(x) 
    cv.glm <- cv.glmnet(x = x, y = y, family = "binomial", alpha = glmnet.alpha, type.measure = "deviance")
    cl <- makeCluster(4)
    stabsel <- function(i) {
        cat("+")
        b_sort <- sort(sample(1 : n, round(3 * n/4)))
        out <- glmnet(x[b_sort, ], y[b_sort], family = "binomial",
                      lambda=cv.glm$lambda.1se, alpha = glmnet.alpha, standardize = FALSE)
        return(tabulate(which(out$beta[, ncol(out$beta)] != 0), p))
    }
    clusterEvalQ(cl, expr = c(library(glmnet)))
    clusterExport(cl, c('stabsel', 'num.runs', 'n', 'glmnet', 'cv.glm', 'x', 'y', 'p'),
                  envir = environment())
    res.cum <- Reduce("+", parLapply(cl, 1 : num.runs, stabsel))
    stopCluster(cl)
    prob.sel <- res.cum / num.runs
    feature.stabsel <- features[prob.sel >= thres]
    return(feature.stabsel)
}

fitModel <- function(train.x, train.y, mdl.name) {
    if (mdl.name == "lasso" || mdl.name == "elasticnet") {
        a <- ifelse(mdl.name == "lasso", yes = 1, no = 0.5) # if elastic-net, then alpha = 0.5
        feature.sel <- extractStableFeatures(x = train.x, y = train.y, thres = 0.5, num.runs = 2e3, glmnet.alpha = a)
        train.x <- train.x[, feature.sel]
        mdl.fit <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = a, type.measure = "deviance")
    } else if (mdl.name == "ridge") {
        mdl.fit <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = 0, type.measure = "deviance")
    } else if (mdl.name == "randomforest") {
        mdl.fit <- randomForest(x = train.x, y = train.y)
    } else {
        stop(paste0("unknown model name to fit: ", mdl.name, " (acceptable options: lasso, elasticnet, ridge, randomforest)"))
    }
    return(mdl.fit)
}
