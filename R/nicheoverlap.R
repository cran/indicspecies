nicheoverlap <- function (P1, P2 = NULL, D = NULL, q = NULL, mode = "multiple", Np1 = NULL, 
            Np2 = NULL, Nq = NULL, nboot = 1000, alpha = 0.05) {
  MODES <- c("single", "multiple", "pairwise")
  mode <- match.arg(mode, MODES)
  if(mode =="multiple" || mode =="single") {
    if(is.null(P2)) stop("P2 cannot be null in mode 'multiple' or 'single'")
    if (!inherits(P1, "data.frame") || !inherits(P2, "data.frame")) stop("P1 and P2 should be dataframes")
    if (mode == "multiple" && (nrow(P1) != nrow(P2))) stop("Resource use dataframes do not have the same number of rows")
    if (ncol(P1) != ncol(P2)) stop("Resource use dataframes do not have the same number of columns (resources)")  
    if (!is.null(Np1) && !is.null(Np2)) {
      if (!inherits(Np1, "numeric") || !inherits(Np2, "numeric")) stop("Np1 and Np2 should be numeric vectors")
      Np1 = as.integer(Np1)
      Np2 = as.integer(Np2)
    }
    if (!is.null(D)) {
      if (!inherits(D, "dist")) 
        stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      if (ncol(P2) != nrow(D)) stop("The number of columns in P2 must be equal to the number of items in D")
      D <- as.dist(D)
    }
  } else if(mode=="pairwise") {
    if (!inherits(P1, "data.frame")) stop("P1 should be a dataframe")
    if (!is.null(D)) {
      if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      D <- as.dist(D)
    }    
  }
  #Check 'q' the availability of resources
  if (!is.null(q)) {
    if (length(q) != ncol(P1)) stop("The number of items in q must be equal to the number of columns in P1 and P2")
    q = q/sum(q)
  } else {
    q = rep(1/ncol(P1), ncol(P1))
  }
  #If no distance matrix was supplied, generate one (equidistant resources)
  if (is.null(D)) D <- as.dist((matrix(1, ncol(P1), ncol(P1)) - diag(rep(1, ncol(P1)))))
  
  #Internal functions
  nichevar1 <- function(f, D) {
      if (is.na(sum(f))) v <- NA
      else if (sum(f) < 1e-16) v <- 0
      else v <- (f %*% (as.matrix(D)^2) %*% f)/(2 * (sum(f)^2))
      return(v)
  }
  getF <- function(p, q = NULL) {
    if (!is.null(q)) { 
      a = p/q
      return(a/sum(a))
    } else {
      return(p/sum(p))
    }
  }
  overlap1 <- function(f1, f2, D) {
    if (is.na(sum(f1)) || is.na(sum(f2))) o <- NA
    else if (sum(f1) < 1e-16 || sum(f2) < 1e-16) o <- 0
    else {
      o <- (1 - ((f1 %*% (as.matrix(D)^2) %*% f2)/(sum(f2) * sum(f1))))/sqrt((1 - 2 * nichevar1(f1, D)) * (1 - 2 * nichevar1(f2, D)))
    }
    return(o)
  }
  
    
  #Calculate overlap between each row of P1 and the corresponding row in P2
  if (mode == "multiple") {
    if (!is.null(Np1) && !is.null(Np2)) nc = 3 #If we have to calculate confidence intervals
    else nc = 1
    O <- as.data.frame(matrix(0, nrow = nrow(P1), ncol = nc))
    for (i in 1:nrow(P1)) rownames(O)[i] <- paste(row.names(P1)[i], "vs", row.names(P2)[i])
    for (i in 1:nrow(P1)) {
      pi1 = as.numeric(P1[i, ])
      pi2 = as.numeric(P2[i, ])
      O[i, 1] <- overlap1(getF(pi1, q), getF(pi2, q), D)
      if (!is.null(Np1) && !is.null(Np2)) {
        BO = vector("numeric", length = nboot)
        O[i, 2] = NA
        O[i, 3] = NA
        if (!is.na(sum(pi1)) && !is.na(sum(pi2))) {
          bsamp1 = rmultinom(nboot, Np1[i], getF(pi1))
          bsamp2 = rmultinom(nboot, Np2[i], getF(pi2))
          if (!is.null(Nq)) qsamp = rmultinom(nboot, Nq, q)
          for (b in 1:nboot) {
            if (!is.null(Nq)) BO[b] = overlap1(getF(bsamp1[, b], qsamp[, b]), getF(bsamp2[, b], qsamp[, b]), D)
            else BO[b] = overlap1(getF(bsamp1[, b], q), getF(bsamp2[, b], q), D)
          }
          #Some NA may appear because of zeroes in qsamp
          BO = BO[!is.na(BO)]
          #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
          z0 = qnorm(sum(BO < O[i, 1])/length(BO))
          lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
          uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
          if (lj > 0 && uj > 0 && lj != uj) {
            sbo = sort(BO)
            O[i, 2] = sbo[lj]
            O[i, 3] = sbo[uj]
          }
        }
      }
    }
    if (nc == 1) names(O) <- "O"
    else names(O) <- c("O", "LC", "UC")
    return(O)
  }
  
  #P1 are different resource use assessments of a single entity (the same for P2)
  if (mode == "single") {
    O <- as.data.frame(matrix(NA, nrow = 1, ncol = 3))
    rownames(O) <- "Overlap"
    O[1, 1] <- overlap1(getF(colSums(P1, na.rm=TRUE), q), getF(colSums(P2, na.rm=TRUE), q), D)
    BO = vector("numeric", length = nboot)
    if (!is.null(Nq)) qsamp = rmultinom(nboot, Nq, q)
    for (b in 1:nboot) {
      p1samp = colSums(P1[sample(1:nrow(P1), replace = TRUE), ], na.rm=TRUE)
      p2samp = colSums(P2[sample(1:nrow(P2), replace = TRUE), ], na.rm=TRUE)
      if (!is.null(Nq)) BO[b] = overlap1(getF(p1samp, qsamp[, b]), getF(p2samp, qsamp[, b]), D)
      else BO[b] = overlap1(getF(p1samp, q), getF(p2samp, q), D)
    }
    BO = BO[!is.na(BO)]
    z0 = qnorm(sum(BO < O[1, 1])/length(BO))
    lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
    uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
    if (lj > 0 && uj > 0 && lj != uj) {
      sbo = sort(BO)
      O[1, 2] = sbo[lj]
      O[1, 3] = sbo[uj]
    }
    names(O) <- c("O", "LC", "UC")
    return(O)
  }
  
  #Only P1 is used. Calculate overlap between pairs of rows in P1. No confidence intervals are calculated
  if(mode=="pairwise") {
    O <- matrix(1, nrow = nrow(P1), ncol = nrow(P1))
    rownames(O)<-rownames(P1)
    colnames(O)<-rownames(P1)
    if(!is.null(Np1)) {
      LC <- O
      UC <- O
    }
    for (i in 1:(nrow(P1)-1)) {
      for (j in (i+1):nrow(P1)) {
        pi = as.numeric(P1[i, ])
        pj = as.numeric(P1[j, ])
        O[i,j] <- overlap1(getF(pi, q), getF(pj, q), D)
        O[j,i] <- O[i,j]
        if (!is.null(Np1)) {
          BO = vector("numeric", length = nboot)
          if (!is.na(sum(pi)) && !is.na(sum(pj))) {
            bsampi = rmultinom(nboot, Np1[i], getF(pi))
            bsampj = rmultinom(nboot, Np1[j], getF(pj))
            if (!is.null(Nq)) qsamp = rmultinom(nboot, Nq, q)
            for (b in 1:nboot) {
              if (!is.null(Nq)) BO[b] = overlap1(getF(bsampi[, b], qsamp[, b]), getF(bsampj[, b], qsamp[, b]), D)
              else BO[b] = overlap1(getF(bsampi[, b], q), getF(bsampj[, b], q), D)
            }
            BO = BO[!is.na(BO)]
            z0 = qnorm(sum(BO < O[i,j])/length(BO))
            lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
            uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
            if (lj > 0 && uj > 0 && lj != uj) {
              sbo = sort(BO)
              LC[i, j] = LC[j, i] =sbo[lj]
              UC[i, j] = UC[j, i] =sbo[uj]
            }
          }
        }
      }
    }
    if(!is.null(Np1)) {
      return(list(O=O, LC=LC, UC = UC))
    }
    else return(O)
  }
}

