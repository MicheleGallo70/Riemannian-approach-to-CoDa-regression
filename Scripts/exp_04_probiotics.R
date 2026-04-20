# Compositional predictor and Scalar  response
df<-data.frame(probiotics)
# Extract the compositional response (Y) and scalar predictor (x_raw)
# to match the function arguments in the rest of the script
X <- as.matrix(df[, 1:130])

# Ensure compositions sum to 1 (closure operation), just to be safe
X <- X / rowSums(X) 

# Natural log of cholesteryl oleate is commonly used for this dataset as per standard CoDa literature

y_raw <- df$y_raw  

res<-ipf_coda_reg(
  predictor=X,
  response=y_raw,
  type   = "comp_pred")