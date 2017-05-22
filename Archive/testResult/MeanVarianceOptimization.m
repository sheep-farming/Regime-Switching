function[frontWts, frontRet, frontVol] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method)

    numAssets = length(ExpRet);
    V0 = zeros(1, numAssets);
    V1 = ones(1, numAssets);

    % Set constraints
    
    UB = ones(numAssets,1);
    
    switch method
        case 'unconstr'
            LB = ones(numAssets,1)*(-10);
            
        case 'constr'
            LB = zeros(numAssets,1);
    end
    
    % Define inequality constraints
    
    A1 = zeros(1,numAssets); 
    A2 = zeros(1,numAssets); 
    A = [A1; A2];
    b = [1; 1];
    
    
    % Find minimum variance return
    minVarWts = quadprog(CovMat, V0, A, b, V1, 1, LB);
    minVarRet = minVarWts'*ExpRet;
    minVarVol = sqrt(minVarWts'*CovMat*minVarWts);
    
    % Find maximum return
    maxRetWts = quadprog([], -ExpRet, A, b, V1, 1, LB);
    maxRet = maxRetWts'*ExpRet;
    
    
    % Calculate frontier portfolios
    targetRet = linspace(minVarRet, maxRet, numPort);
    frontRet = zeros(numPort,1);
    frontVol = zeros(numPort,1);
    frontWts = zeros(numAssets, numPort);
    
    frontRet(1,1) = minVarRet;
    frontVol(1,1) = minVarVol;
    frontWts(:,1) = minVarWts;
    
    Aeq = [V1; ExpRet'];

    
    for i = 2:numPort
        beq = [1; targetRet(i)];
        weights = quadprog(CovMat, V0, A, b, Aeq, beq, LB, UB);
        
        
        frontRet(i,1) = weights'*ExpRet;
        frontVol(i,1) = sqrt(weights'*CovMat*weights);
        frontWts(:,i) = weights;
    end
   
    




