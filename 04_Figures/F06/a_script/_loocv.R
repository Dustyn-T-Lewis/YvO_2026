# Nested LOOCV logistic classifiers on the k features most correlated with the
# label in each training fold. Shared by _supp_prepare_roc.R and
# _supp_multivariate.R; features are re-ranked inside every fold, so selection
# never sees the held-out subject.

rank_features <- function(x, y) {
  r <- abs(cor(x, y))
  r[is.na(r)] <- 0
  names(sort(r[, 1], decreasing = TRUE))
}

fit_predict <- function(train_x, train_y, test_x, feats) {
  fit <- tryCatch(
    suppressWarnings(glm(y ~ .,
      family = binomial,
      data = cbind(y = train_y, as.data.frame(train_x[, feats, drop = FALSE]))
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  predict(fit, type = "response", newdata = as.data.frame(test_x[, feats, drop = FALSE]))
}

# k is chosen per outer fold by inner-LOOCV deviance.
run_topk_loocv <- function(labels, x, k_range = 2:5) {
  x <- as.matrix(x)
  n <- length(labels)
  probs <- numeric(n)
  ks <- integer(n)
  selected <- vector("list", n)

  for (i in seq_len(n)) {
    tr_x <- x[-i, , drop = FALSE]
    tr_y <- labels[-i]
    ranked <- rank_features(tr_x, tr_y)

    inner <- setNames(numeric(length(k_range)), k_range)
    for (k in k_range) {
      feats <- ranked[seq_len(min(k, length(ranked)))]
      dev <- 0
      for (j in seq_along(tr_y)) {
        p <- fit_predict(tr_x[-j, , drop = FALSE], tr_y[-j], tr_x[j, , drop = FALSE], feats)
        if (is.null(p)) {
          dev <- dev + log(2)
          next
        }
        p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
        dev <- dev - (tr_y[j] * log(p) + (1 - tr_y[j]) * log(1 - p))
      }
      inner[as.character(k)] <- dev
    }

    ks[i] <- as.integer(names(which.min(inner)))
    selected[[i]] <- ranked[seq_len(ks[i])]
    probs[i] <- fit_predict(tr_x, tr_y, x[i, , drop = FALSE], selected[[i]]) %||% 0.5
  }
  list(probs = probs, ks = ks, selected = selected)
}

# Fixed k, no inner loop: fast enough to run under label permutation.
fast_loocv_auc <- function(labels, x, k) {
  x <- as.matrix(x)
  probs <- vapply(seq_along(labels), function(i) {
    tr_x <- x[-i, , drop = FALSE]
    feats <- rank_features(tr_x, labels[-i])[seq_len(min(k, ncol(x)))]
    fit_predict(tr_x, labels[-i], x[i, , drop = FALSE], feats) %||% 0.5
  }, numeric(1))
  tryCatch(as.numeric(pROC::auc(pROC::roc(labels, probs, quiet = TRUE))),
    error = function(e) 0.5
  )
}
