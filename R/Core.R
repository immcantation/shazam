#### Transformation functions ####

# Converts a matrix to a vector
#
# \code{clearConsole} clears the console.
# 
# @examples
# # Generate a sample mutations_array 
# sample_matrix <- matrix(sample(20,4),nrow=2, dimnames=list( c("cdr","fwr"), c("r","s") ))
# collapseMatrixToVector(sample_matrix)
#
collapseMatrixToVector <- function(mat, byrow = FALSE){    
    # Get the row and column names
    rnames <- rownames(mat)
    cnames <- colnames(mat)
    if (is.null(rnames)) { rnames <- paste0("Row", 1:nrow(mat)) }
    if (is.null(cnames)) { cnames <- paste0("Column", 1:ncol(mat)) }
    # Combine the row and columns names
    combinedNames <- outer(rnames, cnames, paste, sep = "_")
    
    # Collapse the matrix to a vector
    if (byrow) {
        collapsed_mat <- c(t(mat))
        names(collapsed_mat) <- c(t(combinedNames))
    }
    else{
        collapsed_mat <- c(mat)
        names(collapsed_mat) <- c(combinedNames)
    }
    return(collapsed_mat)
}

# Convert columns to uppercase
#
# @param   data     data.frame to modify.
# @param   columns  vector of column names to transform to uppercase.
# @return  The input data.frame with all entries in \code{columns} transformed 
#          to uppercase.
toupperColumns <- function(data, columns) {
    data <- mutate_at(data, columns, toupper)
    return(data)
}
