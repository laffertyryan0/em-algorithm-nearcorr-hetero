function out = ensureValidNearCorrInput(inputMatrix,tol)
    % This returns a matrix that is made from the input but
    % is ready to be ingested by the nearcorr function
    % This means clipping all entries to be in the range [-1,1]
    % and ensuring the matrix is symmetric with diagonal 1.
    % We hope these violations are small, and if not we should 
    % trigger a warning. 
    matrix = inputMatrix;
    k = length(matrix);
    matrix = matrix - diag(diag(matrix)) + eye(k);
    matrix((matrix >= 1) & ~eye(k)) = 1; 
    matrix((matrix <= -1) & ~eye(k)) = -1;
    out = matrix;
    if max(abs(out-inputMatrix),[],'all')>tol
        warning(...
            "inputMatrix is not close enough to a valid nearcorr input")
    end
end