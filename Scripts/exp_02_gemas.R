# Scalar predictor and compositional response
df<-as.data.frame(gemas_soils)
# Extract the compositional response (Y) and scalar predictor (x_raw)
# to match the function arguments in the rest of the script
Y <- as.matrix(df[,1:18])

# Ensure compositions sum to 1 (closure operation), just to be safe
Y <- Y / rowSums(Y) 

# Natural log of LOI is commonly used for this dataset as per standard CoDa literature
x_raw <- log(df$LOI)                                   # scalar predictor


res<-ipf_coda_reg(
  predictor=x_raw ,
  response=Y,
  type   = "scalar_pred")