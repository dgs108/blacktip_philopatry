library(vegan)
library(stringr)

ordiR2step_atf <- function (object, scope, Pin = 0.01, R2scope = TRUE, permutations = how(nperm = 499),
                        trace = TRUE, R2permutations = 1000, VIF = 3, ...)
{
  if (is.null(object$terms))
    stop("ordination model must be fitted using formula")
  if (missing(scope))
    stop("needs scope")
  if (inherits(scope, "cca"))
    scope <- delete.response(formula(scope))
  if (!inherits(scope, "formula"))
    scope <- reformulate(scope)
  if (is.null(object$CCA))
    R2.0 <- 0
  else R2.0 <- RsquareAdj(object, permutations = R2permutations,
                          ...)$adj.r.squared
  if (is.list(scope) && length(scope) <= 2L)
    scope <- scope$upper
  if (is.null(scope) || !length(add.scope(object, scope)))
    stop("needs upper 'scope': no terms can be added")
  if (R2scope)
    R2.all <- RsquareAdj(update(object, delete.response(formula(scope))), permutations = R2permutations, ...)
  else R2.all <- list(adj.r.squared = NA)
  if (is.na(R2.all$adj.r.squared) && R2scope)
    stop("the upper scope cannot be fitted (too many terms?)")
  R2.all <- R2.all$adj.r.squared
  anotab <- list()
  R2.previous <- R2.0
  high.vif <- NULL
  repeat {
    if (trace) {
      cat("Step: R2.adj=", R2.previous, "\n")
      cat(pasteCall(formula(object)), "\n")
    }
    #print("object:")
    #print(as.character(object$call)[2])
    #print("scope")
    #print(scope)
    adds <- add.scope(object, scope)
    if (length(adds) == 0)
      break
    R2.adds <- numeric(length(adds))
    adds <- paste("+", adds)
    names(R2.adds) <- adds
    for (trm in seq_along(R2.adds)) {
      fla <- paste(". ~ .", names(R2.adds[trm]))
      R2.tmp <- RsquareAdj(update(object, fla), permutations = R2permutations)$adj.r.squared
      if (!length(R2.tmp))
        R2.tmp <- 0
      R2.adds[trm] <- R2.tmp
    }
    best <- which.max(R2.adds)
    if (trace) {
      out <- sort(c(`<All variables>` = R2.all, `<none>` = R2.previous,
                    R2.adds), decreasing = TRUE)
      out <- as.matrix(out)
      colnames(out) <- "R2.adjusted"
      print(out)
      cat("\n")
    }
    if (R2.adds[best] > R2.previous && (!R2scope || R2scope &&
                                        R2.adds[best] <= R2.all)) {
      adds <- names(sort(R2.adds, decreasing=T))
      vif.test <- VIF+1
      while(max(vif.test, na.rm=T) > VIF){
        object_tmp <- object
        fla <- paste("~  .", adds[1])
        object_tmp <- update(object_tmp, fla)
        vif.test <- vif.cca(object_tmp)
        if (max(vif.test, na.rm=T) > VIF) {
          if (trace) {
            print(paste("Max VIF =", max(vif.test, na.rm=T)))
            print(paste(str_sub(adds[1], start=3), "removed from consideration"))
          }
          high.vif<- rbind(high.vif, c(str_sub(adds[1], start=3),max(vif.test, na.rm=T),as.character(object$call)[2]))
          scope <- terms.formula(scope)
          if(length(attr(scope,"term.labels")) == 2){break}
          scope <- reformulate(attr(scope, "term.labels")[!attr(scope, "term.labels") %in% str_sub(adds[1], start=3)])
          adds <- adds[-1]
          if(length(adds) == 0){break}
          #print(scope)
        }
      }
      if(length(adds) == 0){break}
      tst <- add1(object, scope = adds[1], test = "permu", permutations = permutations, alpha = Pin, trace = FALSE)
      if (trace) {
        print(tst[-1, ])
        cat("\n")
      }
      if (is.na(tst[, "Pr(>F)"][2])){
        print(paste(adds[1], "has no constrained component"))
        scope <- terms.formula(scope)
        if(length(attr(scope,"term.labels")) == 2){break}
        scope <- reformulate(attr(scope, "term.labels")[!attr(scope, "term.labels") %in% str_sub(adds[1], start=3)])
        adds <- adds[-1]
      } else if (tst[, "Pr(>F)"][2] <= Pin) {
        fla <- paste("~  .", adds[1])
        object <- update(object, fla)
      } else break
    }
    else {
      break
    }
    R2.previous <- RsquareAdj(object, permutations = R2permutations)$adj.r.squared
    anotab <- rbind(anotab, cbind(R2.adj = R2.previous, tst[2, ]))
  }
  if (NROW(anotab)) {
    if (R2scope)
      anotab <- rbind(anotab, `<All variables>` = c(R2.all, rep(NA, 4)))
    class(anotab) <- c("anova", class(anotab))
    object$anova <- anotab
  }
  if(nrow(high.vif)>0){colnames(high.vif) <- c("Rejected_variable","Max_VIF","Formula")}
  attr(object,'VIF_remove_list') <- high.vif
  object
}