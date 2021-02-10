### Title:    Subroutines for Student-Use Functions


###--------------------------------------------------------------------------###

## Extract the DV name from an lm.fit object
## NOTE: This function only works when lm() is run using the fomula interface.
dvName <- function(x) all.vars(x$terms)[1]

###--------------------------------------------------------------------------###

## Compute the cross-validation error:
getCve <- function(model, data, K, part) {
    ## Loop over K repititions:
    mse <- c()
    for(k in 1 : K) {
        ## Partition data:
        train <- data[part != k, ]
        valid <- data[part == k, ]
        
        ## Fit model, generate predictions, and save the MSE:
        fit    <- lm(model, data = train)
        pred   <- predict(fit, newdata = valid)
        mse[k] <- MSE(y_pred = pred, y_true = valid[ , dvName(fit)])
    }
    ## Return the CVE:
    sum((table(part) / length(part)) * mse)
}

###--------------------------------------------------------------------------###
