# =============================================================================

ipf_coda_reg <- function(
    predictor,
    response,
    type   = c("scalar_pred", "comp_pred"),
    t_end  = 500,
    n_time = 5000,
    tol    = 1e-8,
    eps    = 1e-12
) {
  # =============================================================================
  #  ipf_coda_reg()
  #  Riemannian IPF regression for Compositional Data (CoDa)
  #
  #  type = "scalar_pred"  : scalar predictor x in R,
  #                          compositional response Y in V*(p-1)   [Eq. 28]
  #  type = "comp_pred"    : compositional predictor X in V*(p-1),
  #                          scalar response y in R                [Eq. 31]
  #
  #  Inference is carried out in CLR coordinates:
  #   - scalar_pred : p independent simple regressions (one per CLR component)
  #   - comp_pred   : multiple regression via pseudoinverse of X_clr^T X_clr
  #
  #  Dependencies : deSolve
  # =============================================================================
  # ---- 0. Validate inputs ---------------------------------------------------
  
  if (!requireNamespace("deSolve", quietly = TRUE))
    stop("Package 'deSolve' is required. Install it with install.packages('deSolve').")
  
  .is_tabular    <- function(x) is.matrix(x) || is.data.frame(x)
  .is_scalar_vec <- function(x) {
    if (is.numeric(x) && is.null(dim(x))) return(TRUE)
    if (.is_tabular(x) && ncol(x) == 1L)  return(TRUE)
    FALSE
  }
  
  if (missing(type)) {
    if (.is_tabular(predictor) && .is_scalar_vec(response)) {
      type <- "comp_pred"
      message("Note: 'type' inferred as \"comp_pred\" (tabular predictor, scalar response).")
    } else if (.is_scalar_vec(predictor) && .is_tabular(response)) {
      type <- "scalar_pred"
      message("Note: 'type' inferred as \"scalar_pred\" (scalar predictor, tabular response).")
    }
  }
  type <- match.arg(type)
  
  if (type == "scalar_pred" && .is_tabular(predictor)) {
    hint <- if (.is_scalar_vec(response)) " Did you mean type = \"comp_pred\"?" else ""
    stop("For type = \"scalar_pred\", 'predictor' must be a numeric vector, not a ",
         class(predictor)[1L], " with ", ncol(predictor), " columns.", hint)
  }
  if (type == "comp_pred" && .is_scalar_vec(predictor) && .is_tabular(response)) {
    stop("For type = \"comp_pred\", 'predictor' must be a compositional matrix. ",
         "Did you mean type = \"scalar_pred\"?")
  }
  
  # ---- internal helpers -----------------------------------------------------
  
  closure <- function(v) { v <- pmax(v, eps); v / sum(v) }
  
  clr_vec <- function(v) { lv <- log(pmax(v, eps)); lv - mean(lv) }
  
  clr_mat <- function(M) t(apply(M, 1, clr_vec))
  
  aitchison_mean <- function(M) closure(exp(colMeans(log(pmax(M, eps)))))
  
  center_comp <- function(M) {
    Mc <- sweep(M, 2, aitchison_mean(M), "/")
    Mc / rowSums(Mc)
  }
  
  # =========================================================================
  #  ODE RIGHT-HAND SIDES
  # =========================================================================
  
  .flux_scalar_pred <- function(t, a, params) {
    xc    <- params$xc
    Y_cen <- params$Y_cen
    n     <- nrow(Y_cen)
    p     <- length(a)
    a     <- closure(a)
    ln_a  <- log(a)
    ones  <- rep(1, p)
    acc   <- numeric(p)
    for (i in seq_len(n)) {
      a_xi     <- a^xc[i]
      s_xi     <- sum(a_xi)
      delta_i  <- log(Y_cen[i, ]) - xc[i] * ln_a + log(s_xi)
      acc      <- acc + xc[i] * (a_xi / s_xi - ones) * delta_i
    }
    nabla_a <- acc * a
    list(as.numeric(a * sum(nabla_a) - nabla_a))
  }
  
  .flux_comp_pred <- function(t, a, params) {
    X_clr <- params$X_clr
    y_cen <- params$y_cen
    a     <- closure(a)
    clr_a <- clr_vec(a)
    delta <- y_cen - as.numeric(X_clr %*% clr_a)
    nabla <- as.numeric(t(X_clr) %*% delta)
    proj  <- nabla - mean(nabla)           # (I - 11^T/p) nabla
    a_dot <- a * proj - a * sum(a * proj)
    list(as.numeric(a_dot))
  }
  
  .root_scalar <- function(t, a, p) max(abs(.flux_scalar_pred(t, a, p)[[1]])) - tol
  .root_comp   <- function(t, a, p) max(abs(.flux_comp_pred(t, a, p)[[1]]))   - tol
  
  # =========================================================================
  #  DISPATCH
  # =========================================================================
  
  times_seq <- seq(0, t_end, length.out = n_time)
  
  # -------------------------------------------------------------------------
  if (type == "scalar_pred") {
    
    x <- as.numeric(predictor)
    Y <- as.matrix(response)
    if (length(x) != nrow(Y))
      stop("'predictor' length must equal number of rows of 'response'.")
    
    n <- nrow(Y); p <- ncol(Y)
    x_mean <- mean(x);  xc <- x - x_mean
    y_mu   <- aitchison_mean(Y)
    Y_cen  <- center_comp(Y)
    
    a_init        <- rep(1 / p, p)
    names(a_init) <- colnames(Y)
    
    sol <- deSolve::ode(y = a_init, times = times_seq,
                        func = .flux_scalar_pred,
                        parms = list(xc = xc, Y_cen = Y_cen),
                        rootfun = .root_scalar)
    
    a_final <- closure(as.numeric(sol[nrow(sol), -1]))
    names(a_final) <- colnames(Y)
    
    b_comp <- closure(y_mu / closure(a_final^x_mean))
    names(b_comp) <- colnames(Y)
    
    Y_hat <- t(sapply(x, function(xi) closure(a_final^xi * b_comp)))
    colnames(Y_hat) <- colnames(Y)
    
    # Aitchison residuals
    res_ait <- sapply(seq_len(n), function(i) {
      r <- clr_vec(closure(Y[i, ] / Y_hat[i, ]))
      sqrt(sum(r^2))
    })
    
    # ---- Inference in CLR space --------------------------------------------
    
    clr_a   <- clr_vec(a_final)
    clr_b   <- clr_vec(b_comp)
    clr_Y   <- clr_mat(Y)
    clr_ymu <- clr_vec(y_mu)
    
    df_res  <- n - 1L
    
    # Overall Aitchison R²
    SS_res <- sum(res_ait^2)
    SS_tot <- sum(sapply(seq_len(n), function(i) sum((clr_Y[i, ] - clr_ymu)^2)))
    R2_ait     <- 1 - SS_res / SS_tot
    k          <- p - 1L                  # free slope parameters (CLR constraint)
    R2_adj_ait <- 1 - (1 - R2_ait) * (n - 1) / (n - k - 1)
    
    return(structure(
      list(
        type             = "scalar_pred",
        n                = n, p = p,
        a_slope          = a_final,
        b_intercept      = b_comp,
        clr_a            = clr_a,
        clr_b            = clr_b,
        df_residual      = df_res,
        Y_hat            = Y_hat,
        residuals_ait    = res_ait,
        RMSE_ait         = sqrt(mean(res_ait^2)),
        R2               = R2_ait,
        R2_adj           = R2_adj_ait,
        convergence_time = sol[nrow(sol), "time"],
        ode_path         = sol
      ), class = "ipf_coda_fit"))
  }
  
  # -------------------------------------------------------------------------
  # type == "comp_pred"
  
  X <- as.matrix(predictor)
  y <- as.numeric(response)
  if (nrow(X) != length(y))
    stop("Number of rows of 'predictor' must equal length of 'response'.")
  
  n <- nrow(X); p <- ncol(X)
  y_mean <- mean(y);  y_cen <- y - y_mean
  x_mu   <- aitchison_mean(X)
  X_cen  <- center_comp(X)
  X_clr  <- clr_mat(X_cen)
  
  a_init        <- rep(1 / p, p)
  names(a_init) <- colnames(X)
  
  sol <- deSolve::ode(y = a_init, times = times_seq,
                      func = .flux_comp_pred,
                      parms = list(X_clr = X_clr, y_cen = y_cen),
                      rootfun = .root_comp)
  
  a_final <- closure(as.numeric(sol[nrow(sol), -1]))
  names(a_final) <- colnames(X)
  
  clr_a  <- clr_vec(a_final)
  clr_xm <- clr_vec(x_mu)
  b_int  <- y_mean - sum(clr_a * clr_xm)
  
  X_clr_orig <- clr_mat(X)
  y_hat      <- as.numeric(X_clr_orig %*% clr_a) + b_int
  
  res <- y - y_hat
  SSR <- sum(res^2)
  SST <- sum((y - y_mean)^2)
  R2     <- 1 - SSR / SST
  k      <- p - 1L                                 # free slope parameters
  R2adj  <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
  
  
  df_res    <- n - p
  
  structure(
    list(
      type             = "comp_pred",
      n                = n, p = p,
      a_slope          = a_final,
      b_intercept      = b_int,
      clr_a            = clr_a,
      df_residual      = df_res,
      y_hat            = y_hat,
      residuals        = res,
      RMSE             = sqrt(mean(res^2)),
      R2               = R2,
      R2_adj           = R2adj,
      convergence_time = sol[nrow(sol), "time"],
      ode_path         = sol
    ), class = "ipf_coda_fit")
}