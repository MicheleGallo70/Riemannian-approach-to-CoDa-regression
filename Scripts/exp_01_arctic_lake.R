# Scalar predictor and compositional response
df<-as.data.frame(arctic_lake)
# Extract the compositional response (Y) and scalar predictor (x_raw)
# to match the function arguments in the rest of the script
Y <- as.matrix(df[, c("sand", "silt", "clay")])

# Ensure compositions sum to 1 (closure operation), just to be safe
Y <- Y / rowSums(Y) 

# Natural log of depth is commonly used for this dataset as per standard CoDa literature
x_raw <- log(df$depth)                                   # scalar predictor


res<-ipf_coda_reg(
  predictor=x_raw ,
  response=Y,
  type   = "scalar_pred")