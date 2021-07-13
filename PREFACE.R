#!/usr/bin/env Rscript

version = 'v0.1.2'

# ---
# Functions
# ---

print.help <- function(w = 'all'){
  cat('\nUsage:\n')
  if (w == 'all' | w == 'train'){
    cat('\tRScript PREFACE.R train --config path/to/config.txt --outdir path/to/dir/ [--nfeat (int) --hidden (int) --cpus (int) --femprop --olm --noskewcorrect]\n')
  }
  if (w == 'all' | w == 'predict'){
    cat('\tRScript PREFACE.R predict --infile path/to/infile.bed --model path/to/model.RData [–-json]\n')
  }
  cat('\n')
  quit(save = 'no')
}

unrec.args <- function(w = 'all'){
  cat('\nUnrecongized arguments, have you read the manual at \'https://github.com/CenterForMedicalGeneticsGhent/PREFACE\'?')
  print.help(w)
}

parse.op.arg <- function(args, sub.arg, default){
  if (sub.arg %in% args){
    i = which(args == sub.arg) + 1
    resp = as.integer(args[i])
    args = args[-i]
  } else {
    resp = default
  }
  if (is.na(resp)){
    cat(paste0('Argument \'', sub.arg, '\' requires a value.\n'))
    print.help()
  }
  return(list(resp, args))
}

get.m.diff <- function(v1, v2, abs = T){
  if (!abs){
    return(mean(v1 - v2))
  }
  return(mean(abs(v1 - v2)))
}

get.sd.diff <- function(v1, v2){
  return(sd(abs(v1 - v2)))
}

plot.performance <- function(v1, v2, summary, n.feat, xlab, ylab, path){
  
  png(path, width=7.6, height=2.65, units='in', res=1024)
  par(mar=c(3,1,2,1), mgp=c(1.6, 0.2, 0.2), mfrow=c(1,3), xpd = NA, oma=c(0,3,0,0))
  
  ylim = c(min(summary$importance[2,][summary$importance[2,] != 0]), max(summary$importance[2,]))
  xlim = c(1, ncol(summary$importance))
  
  plot(log(1:ncol(summary$importance)), log(summary$importance[2,]), ylim = log(ylim), xlim = log(xlim), type = 'l', lwd = 2,
       axes = F, ylab = 'Proportion of variance', xlab = 'Principal components', col = color.A, main = 'PCA')
  
  segments(log(n.feat), log(ylim[1]), log(n.feat), log(ylim[2] * .99), lwd = 3, lty = 3, c = color.C)
  text(log(n.feat), log(ylim[2]), 'Number of features', col = color.C, adj = 0.5, cex = 0.8)
  
  label.seq = round(seq(from = xlim[1], to = xlim[2], (xlim[2] - xlim[1])/xlim[2]))
  axis(1, tcl=0.5, at = log(label.seq), labels = label.seq)
  label.seq = seq(from = ylim[1], to = ylim[2], (ylim[2] - ylim[1])/xlim[2])
  axis(2, tcl=0.5, at = log(label.seq), labels = rep('', length(label.seq)), las = 2)
  axis(2, tcl=0.5, at = log(c(label.seq[1], label.seq[length(label.seq)])), labels = c(label.seq[1], label.seq[length(label.seq)]), las = 2)
  
  xlim <- c(0, max(v1))
  ylim <- c(0, max(v2))
  mx <- max(xlim[2], ylim[2])
  
  plot(v1, v2, pch = 16, cex = 0.6, axes = F, xlab = xlab, ylab = ylab,
       xlim = c(0, mx), ylim = c(0, mx),
       main = 'Scatter plot')
  
  legend('topleft', c('OLS fit', 'f(x)=x'), bty = 'n',
         col = c(color.A, color.B), cex = 0.9, text.col = c(color.A, color.B), text.font = 2)
  
  axis(1, tcl=0.5)
  axis(2, tcl=0.5, las = 2)
  
  fit <- coef(lsfit(v1, v2))
  
  par(xpd=F)
  segments(0, 0, mx, mx, lwd = 3, lty = 3, c = color.B)
  segments(0, fit[1], mx, fit[1] + mx * fit[2], lwd = 3, lty = 2, c = color.A)
  par(xpd=NA)
  
  text(0, mx * 1.03,
       paste0('(r = ', signif(cor(v1, v2), 3), ')'),
       cex = 0.9, adj = 0)
  t = hist(v1-v2, max(20,length(v1)/10), axes = F, xlab = paste0(xlab, ' - ', ylab),
           main = 'Histogram', ylab = 'Density',  c = 'black')
  
  mx = max(t$counts)
  segments(0, 0, 0, mx, lwd = 3, lty = 3, c = color.B)
  segments(get.m.diff(v1, v2, abs = F), 0, get.m.diff(v1, v2, abs = F), mx, lwd = 3, lty = 2, c = color.A)
  axis(1, tcl=0.5)
  axis(2, tcl=0.5, las = 2)
  legend('topleft', c('mean error', 'x=0'), bty = 'n',
         col = c(color.A, color.B), cex = 0.9, text.col = c(color.A, color.B), text.font = 2)
  
  text(min(t$breaks), mx * 1.03,
       paste0('(MAE = ', signif(get.m.diff(v1, v2), 3), ' ± ',  signif(get.sd.diff(v1, v2), 3), ')'),
       cex = 0.9, adj = 0)
  
  dev.off()
  return(c(fit[1], fit[2], get.m.diff(v1, v2), get.sd.diff(v1, v2), cor(v1, v2)))
}

train.neural <- function(f, train.nn, hidden){
  tryCatch({
    return(neuralnet(f, train.nn, hidden = hidden, stepmax = 1e6))
  }, warning = function(e) {
    cat(paste0('Neural network did not converge. Re-run and decrease --hidden or optimize --nfeat. Alternatively, use --olm.\n'))
    quit(save = 'no')
  })
}

# ---
# Modules
# ---

train <- function(args){
  
  start.time <- proc.time()
  
  # Additional lib
  
  suppressMessages(library('foreach'))
  suppressMessages(library('doParallel'))
  suppressMessages(library('MASS'))
  suppressMessages(library('irlba'))
  
  # Arg parse
  
  ## Mandatory

  man.args <- c('--config', '--outdir')
  if (length(which(man.args %in% args)) != length(man.args)) unrec.args('train')
  config.file = args[which(args == '--config') + 1]
  out.dir = args[which(args == '--outdir') + 1]
  if(!(file.exists(config.file))){
    cat(paste0('The file \'', config.file, '\' does not exist.\n'))
    quit(save = 'no')
  }
  args = args[args != config.file]
  
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  if(!(file.exists(out.dir))){
    cat(paste0('Could not create directory \'', out.dir, '\'.\n'))
    quit(save = 'no')
  }
  args = args[args != out.dir]
  
  ## Optional
  op.args <- c('--nfeat', '--hidden', '--olm', '--femprop',  '--noskewcorrect', '--cpus')
  
  n.feat <- parse.op.arg(args, '--nfeat', 50)[[1]] ; args <- parse.op.arg(args, '--nfeat', 50)[[2]]
  hidden <- parse.op.arg(args, '--hidden', 2)[[1]] ; args <- parse.op.arg(args, '--hidden', 2)[[2]]
  cpus <- parse.op.arg(args, '--cpus', 1)[[1]] ; args <- parse.op.arg(args, '--cpus', 1)[[2]]
  is.olm = F ; if ('--olm' %in% args) is.olm = T
  skewcorrect = T ; if ('--noskewcorrect' %in% args) skewcorrect = F
  train.gender = c('M') ; if ('--femprop' %in% args) train.gender = c('M', 'F')
  
  ## Others
  
  if(any(!(args %in% c(man.args, op.args)))){
    cat(paste0('Argument(s) \'', paste0(args[!(args %in% c(man.args, op.args))], collapse = '\', \''), '\' not recognized. Will ignore.\n'))
  }
  
  out.dir <- paste0(out.dir, '/')
  
  # Start training
  
  ## Load necessary files
  
  config.file <- read.csv(file = config.file, sep = '\t', header = T,
                          comment.char='', colClasses = c('character', 'character', 'factor', 'numeric'))
  if (length(which(config.file$gender %in% train.gender)) < n.feat){
    cat(paste0('Please provide at least ', n.feat, ' labeled samples.\n'))
    quit(save = 'no')
  }
  
  config.file <- config.file[sample(nrow(config.file)),]
  
  training.frame <- read.table(config.file$filepath[1], header = T, sep = '\t')
  training.frame <- training.frame[,colnames(training.frame) %in% c('chr', 'start', 'end')]
  training.frame <- training.frame
  
  registerDoParallel(cpus)
  training.frame.sub <- foreach(i = 1:nrow(config.file), .combine = 'cbind') %dopar% {
    sample <- config.file$ID[i]
    cat(paste0('Loading sample ', sample, ' | ', nrow(config.file) - i, '/', nrow(config.file), ' remaining ...\n'))
    bin.table <- fread(config.file$filepath[i], header = T, sep = '\t')
    return(suppressWarnings(as.numeric(bin.table$ratio[bin.table$chr != 'Y'])))
  }

  colnames(training.frame.sub) <- config.file$ID
  X.ratios <- as.data.frame(training.frame.sub['X' == training.frame$chr, ])
  X.ratios <- 2 ** colMeans(X.ratios, na.rm = T)

  cat(paste0('Creating training frame ...\n'))
  
  training.frame <- cbind(training.frame[!(training.frame$chr %in% exclude.chrs),],
                          training.frame.sub[!(training.frame$chr %in% exclude.chrs),])
  training.frame.t <- t(training.frame[4:ncol(training.frame)])
  colnames(training.frame.t) <- paste0(training.frame$chr, ':', training.frame$start, '-', training.frame$end)
  
  training.frame <-  as.data.frame(training.frame.t)
  rm(training.frame.t)
  
  training.frame <- training.frame[,colSums(is.na(training.frame)) < nrow(config.file) * 0.01]
  possible.features <- colnames(training.frame)
  mean.features <- colMeans(training.frame, na.rm = T)
  
  na.index <- which(is.na(training.frame), arr.ind=TRUE)
  if (length(na.index[,2])) training.frame[na.index] <- mean.features[na.index[,2]]
  
  cat(paste0('Remaining training features after \'NA\' filtering: ', length(possible.features), '\n'))
  
  ## Predictive modeling
  
  dir.create(paste0(out.dir, 'training_repeats'), showWarnings = FALSE, recursive = TRUE)
  
  repeats = 10
  test.percentage = 1/repeats
  
  test.number = length(which(config.file$gender %in% train.gender)) * test.percentage
  
  max.feat <- length(which(config.file$gender %in% train.gender)) - as.integer(test.number) - 1
  if (n.feat > max.feat){
    cat(paste0('Too few samples were provided for --nfeat ', n.feat, ', using --nfeat ', max.feat, '\n'))
    n.feat <- max.feat
  }
  
  oper <- foreach(i = 1:repeats) %dopar% {

    cat(paste0('Model training | Repeat ', i,'/', repeats, ' ...\n'))
    
    test.index.overall <- which(config.file$gender %in% train.gender)[(as.integer((i-1)*test.number) + 1):as.integer(i*test.number)]
    train.index.overall <- sort(which(config.file$gender %in% train.gender)[!((which(config.file$gender %in% train.gender)) %in% test.index.overall)])
    train.index.subset <- sort(which(config.file$gender[-test.index.overall] %in% train.gender))
    
    cat(paste0('\tExecuting principal component analysis ...\n'))
    pca.train <- suppressWarnings(prcomp_irlba(training.frame[-test.index.overall,],
                              n = min(n.feat * 10, nrow(training.frame[-test.index.overall,]) - 1), scale. = F))
    
    X.train <- as.matrix(pca.train$x[train.index.subset, ])
    Y.train <- as.matrix(config.file$FF[train.index.overall], ncol = 1)
    X.test <- as.matrix(scale(training.frame[test.index.overall,], pca.train$center, pca.train$scale) %*% pca.train$rotation)
    Y.test <- as.matrix(config.file$FF[test.index.overall], ncol = 1)
    
    if (is.olm){
      cat(paste0('\tTraining ordinary linear model ...\n'))
      model <- glmnet(x = X.train[,1:n.feat], y = Y.train, family='gaussian', lambda = 0)
      prediction = as.numeric(predict.glmnet(model, X.test[,1:n.feat]))
    } else {
      train.nn <- X.train[,1:n.feat]
      train.nn <- cbind(train.nn, Y.train)
      colnames(train.nn) <- c(colnames(train.nn)[1:(ncol(train.nn) - 1)], 'FF')
      f <- paste(colnames(X.train)[1:(ncol(train.nn) - 1)], collapse=' + ')
      f <- paste('FF ~',f)
      f <- as.formula(f)
      cat(paste0('\tTraining neural network ...\n'))
      model <- train.neural(f, train.nn, hidden)
      prediction = as.numeric(compute(model, X.test[,1:n.feat])$net.result)
    }
    
    info <- plot.performance(prediction, Y.test, summary(pca.train), n.feat, 'PREFACE (%)', 'FF (%)', paste0(out.dir, 'training_repeats/', 'repeat_', i,'.png'))
    
    results <- list()
    results$intercept <- as.numeric(info[1])
    results$slope <- as.numeric(info[2])
    results$prediction <- prediction
    return(results)
  }
  stopImplicitCluster()

  predictions <- c()
  for(rep in 1:repeats){
    predictions <- c(predictions, oper[[rep]]$prediction)
  }
  
  if (skewcorrect){
    p <- sample(length(predictions))[1:(length(predictions)/4)]
    fit <- coef(lsfit(predictions[p], config.file$FF[config.file$gender %in% train.gender][p]))
    the.intercept <- fit[1]
    the.slope <- fit[2]
  } else {
    the.intercept <- 0
    the.slope <- 1
  }
  
  ## FFX
  
  png(paste0(out.dir, 'FFX.png'), width=5, height=2.8, units='in', res=1024)
  par(mar=c(2.7,2,0,0.3), mgp=c(1.5, 0.2, 0.2), mfrow=c(1,2), xpd = NA, oma=c(0,1.5,0,0))
  
  v1 <- config.file$FF[config.file$gender == 'M']
  v2 <- X.ratios[config.file$gender == 'M']
  plot(v1, v2, pch = 16, cex = 0.4, axes = F, xlab = 'FF (%)', ylab = 'μ(ratio X)',
       xlim = c(0, max(v1)))
  fit <- rlm(v2 ~ v1)
  
  r2wls <- function(x){
    SSe <- sum(x$w*(x$resid)^2)
    observed <- x$resid+x$fitted
    SSt <- sum(x$w*(observed-weighted.mean(observed,x$w))^2)
    value <- 1-SSe/SSt;
    return(value);
  }
  
  r <- r2wls(fit) ** 0.5 ; fit <- coef(fit)
  
  par(xpd = F)
  segments(min(v1), fit[1] + min(v1) * fit[2], max(v1), fit[1] + max(v1) * fit[2], lwd = 3, lty = 2, c = color.A)
  par(xpd = NA)
  
  axis(1, tcl=0.5)
  axis(2, tcl=0.5, las = 2)
  
  legend('topright', legend = c('RLS fit', 'f(x)=x', paste0('(wr = ', signif(r, 4), ')')),
         bty = 'n', lty = c(2, 3, -1), col = c(color.A, color.B, 'black'), cex = 0.7, text.col = c(color.A, color.B, 'black'),
         text.font = c(2, 2, 1))
  
  
  v2 <- (v2 - fit[1]) / fit[2]
  plot(v1, v2, pch = 16, cex = 0.4, axes = F, xlab = 'FF (%)', ylab = 'FFX (%)',
       xlim = c(0, max(v1)))
  segments(min(v1), min(v1), max(v1), max(v1), lwd = 3, lty = 3, c = color.B)
  
  
  axis(1, tcl=0.5)
  axis(2, tcl=0.5, las = 2)
  
  dev.off()
  
  the.intercept.X <- fit[1]
  the.slope.X <- fit[2]
  
  ## Output final model & accuracy statistics
  
  predictions <- the.intercept + the.slope * predictions
  
  cat(paste0('Executing final principal component analysis ...\n'))
  pca.train <- suppressWarnings(prcomp_irlba(training.frame, n = min(n.feat * 10, nrow(training.frame) - 1), scale. = F))
  X.train <- as.matrix(pca.train$x[which(config.file$gender %in% train.gender), ])
  Y.train <- as.matrix(config.file$FF[which(config.file$gender %in% train.gender)], ncol = 1)
  
  if (is.olm){
    cat(paste0('Training final ordinary linear model ...\n'))
    model <- glmnet(x = X.train[,1:n.feat], y = Y.train, family='gaussian', lambda = 0)
  } else {
    train.nn <- X.train[,1:n.feat]
    train.nn <- cbind(train.nn, Y.train)
    colnames(train.nn) <- c(colnames(train.nn)[1:(ncol(train.nn) - 1)], 'FF')
    f <- paste(colnames(X.train)[1:(ncol(train.nn) - 1)], collapse=' + ')
    f <- paste('FF ~',f)
    f <- as.formula(f)
    cat(paste0('Training final neural network ...\n'))
    model <- train.neural(f, train.nn, hidden)
  }
  
  train.config.file <- config.file[config.file$gender %in% train.gender, ]
  
  info <- plot.performance(predictions, train.config.file$FF, summary(pca.train),
                           n.feat, 'PREFACE (%)', 'FF (%)', paste0(out.dir, 'overall_performance.png'))
  
  index.10 <- which(train.config.file$FF < 20 - predictions)
  deviations.10 <- abs(predictions[index.10] - train.config.file$FF[index.10])
  
  deviations <- abs(predictions - train.config.file$FF)
  outliers.index <- which(deviations > info[4] + 3 * info[5])
  outliers <- train.config.file$ID[outliers.index]
  outlier.values <- deviations[outliers.index]
  outliers.values.noabs <- c(predictions - train.config.file$FF)[outliers.index]
  
  sink(paste0(out.dir, 'training_statistics.txt'))
  cat('PREFACE - PREdict FetAl ComponEnt\n\n')
  if (length(outliers) != 0){
    cat(paste0('Below, some of the top candidates for outlier removal are listed.\n',
      'If you know some of these are low quality/have sex aberrations (when using FFY as response variable), remove them from the config file and re-run.\n',
      'Avoid removing other cases, as this will result in inaccurate performance statistics and possible overfitting towards irrelevant models.\n\n'))
    cat('_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_\n')
    cat('ID\tFF (%) - PREFACE (%)\n')
    for (i in rev(order(outlier.values))){
      cat(paste0(outliers[i], '\t', outliers.values.noabs[i], '\n'))
    }
    cat('_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n')
  }
  elapsed.time <- proc.time() - start.time
  cat(paste0('Training time: ', as.character(round(elapsed.time[3])), ' seconds\n'))
  cat(paste0('Overall correlation (r): ', info[5], '\n'))
  cat(paste0('Overall mean absolute error (MAE): ', info[3], ' ± ', info[4], '\n'))
  cat(paste0('FF < 10% mean absolute error (MAE): ', mean(deviations.10), ' ± ', sd(deviations.10), '\n'))
  cat(paste0('Correction for skew: \n'))
  cat(paste0('\tIntercept: ', the.intercept, '\n'))
  cat(paste0('\tSlope: ', the.slope, '\n\n'))
  cat(paste0('Do not forget to verify whether the \'--nfeat\' parameter captures the first \'random\' phase (and not too much of the \'non-random\' phase) at \'', out.dir, 'overall_performance.png\'.\n',
             'If you believe this parameter is not located in an optimal position, decrease/increase \'--nfeat\' and re-run.'))
  sink()
  
  pca.center <- pca.train$center
  pca.scale <- pca.train$scale
  pca.rotation <- pca.train$rotation[,1:n.feat]
  
  if(is.olm){
    model <- model[which(names(model) %in% c('beta', 'lambda', 'a0', 'offset'))]
  } else {
    model <- model[which(names(model) %in% c('linear.output', 'weights', 'model.list'))]
  }

  save('n.feat', 'mean.features', 'possible.features', 'pca.center',
       'pca.scale', 'pca.rotation', 'model', 'the.intercept',
       'the.slope', 'the.intercept.X', 'the.slope.X', 'is.olm',
       file = paste0(out.dir, 'model.RData'))
  
  cat(paste0('Finished! Consult \'', out.dir, 'training_statistics.txt\' to analyse your model\'s performance.\n'))
}

predict <- function(args){
  
  # Arg parse
  
  man.args <- c('--infile', '--model')
  if (length(which(man.args %in% args)) != length(man.args)) unrec.args('predict')
  in.file = args[which(args == '--infile') + 1]
  model = args[which(args == '--model') + 1]
  if(!(file.exists(in.file))){
    cat(paste0('The file \'', in.file, '\' does not exist.\n'))
    quit(save = 'no')
  }
  if(!(file.exists(model))){
    cat(paste0('The file \'', model, '\' does not exist.\n'))
    quit(save = 'no')
  }
  args <- args[!(args %in% c(in.file, model))]
  
  op.args <- c('--json')
  if ('--json' %in% args){
    json = as.character(args[which(args == '--json') + 1])
    if (!is.na(json)){
      args = args[args != json]
    }
  } else {
    json = ''
  }
  
  if(any(!(args %in% c(man.args, op.args, in.file, model)))){
    cat(paste0('Argument(s) \'', paste0(args[!(args %in% c(man.args, op.args))], collapse = '\', \''), '\' not recognized. Will ignore.\n'))
  }
  
  # Predict

  load(model)
  model$act.fct <- function (x) {1/(1 + exp(-x))}
  bin.table <- fread(in.file, header = T, sep = '\t')
  feat.id <- paste0(bin.table$chr, ':', bin.table$start, '-', bin.table$end)
  X.ratio <- 2 ** mean(as.numeric(bin.table$ratio['X' == bin.table$chr]), na.rm = T)
  ratio <- as.numeric(bin.table$ratio[which(feat.id %in% possible.features)])
  
  FFX <- (X.ratio - the.intercept.X) / the.slope.X
  
  if (any(is.na(ratio))){
    ratio[is.na(ratio)] <- mean.features[is.na(ratio)]
  }
  
  projected.ratio <- as.matrix(scale(t(ratio), pca.center, pca.scale) %*% pca.rotation)
  if (is.olm){
    prediction <- as.numeric(predict.glmnet(model, projected.ratio))
  } else{
    prediction = as.numeric(compute(model, projected.ratio)$net.result)
  }
  
  prediction <- the.intercept + the.slope * prediction
  
  json.dict <- paste0('{\"FFX\": ', FFX / 100, ', \"PREFACE\": ', prediction / 100, '}')
  
  if (is.na(json)){
    cat(json.dict)
    cat('\n')
  } else {
    if (json == ''){
    cat(paste0('FFX = ', signif(FFX, 4), '%\n'))
    cat(paste0('PREFACE = ', signif(prediction, 4), '%\n'))
    } else {
    sink(paste0(json))
    cat(json.dict)
    sink()
    }
  }
}

# ---
# overall lib
# ---

suppressMessages(library('data.table'))
suppressMessages(library('glmnet'))
suppressMessages(library('neuralnet'))

# ---
# param
# ---

exclude.chrs <- c('13', '18', '21', 'X', 'Y')

color.A = rgb(141, 209, 198, maxColorValue = 255)
color.B = rgb(227, 200, 138, maxColorValue = 255)
color.C = rgb(200, 120, 120, maxColorValue = 255)

# ---
# Main
# ---

set.seed(1)
args <- commandArgs(trailingOnly = TRUE)
if ('--help' %in% args | '--h' %in% args) print.help()
if ('--version' %in% args | '--v' %in% args){cat(paste0(version), '\n') ; quit(save = 'no')}
if (length(args) == 0) unrec.args()
if (!(args[1] %in% c('train', 'predict'))){
  unrec.args()
}
if (args[1] == 'train') train(args[-1])

if (args[1] == 'predict') predict(args[-1])
