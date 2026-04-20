# Compositional predictor and Scalar  response
df<-data.frame(coxite)
# Extract the compositional response (Y) and scalar predictor (x_raw)
# to match the function arguments in the rest of the script
X <- as.matrix(df[, 1:6])

# Ensure compositions sum to 1 (closure operation), just to be safe
X <- X / rowSums(X) 

# log(df$porosity)/ (100 - df$porosity) is commonly used for this dataset as per standard CoDa literature

y_raw <- log(df$porosity)/ (100 - df$porosity)     

res<-ipf_coda_reg(
  predictor=X,
  response=y_raw,
  type   = "comp_pred")