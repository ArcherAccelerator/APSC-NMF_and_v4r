clear; clc;
%% NMF Source Apportionment Implementation in MATLAB
% This script performs Non-Negative Matrix Factorization (NMF) for source 
% apportionment of heavy metal concentrations in PM2.5 samples.
% Key features: Sequential planning constraints on factor loadings, data 
% normalization, and comparison with PMF results.
% Author: Yuqinqin (Modified for clarity and professionalism, 2025)
% Date: November 2, 2024 (Original); December 16, 2025 (Revised)

%% Configuration Parameters
numSamples = 95;          % Number of environmental samples
numMetals = 13;           % Number of heavy metals
numFactors = 4;           % Number of source factors (adjustable)
maxIterations = 1000;     % Maximum NMF iterations
numReplicates = 100;      % Number of NMF replicates for robustness

% File paths (parameterize for portability)
inputDataPath = '..\Case_data\SourceApportionmentData.xls';
pmfFactorScoresPath = '..\Case_data\PMF_G_Matrix.xlsx';
pmfFactorLoadingsPath = '..\Case_data\PMF_F_Matrix.xlsx';
sequenceConstraintsPath = '..\Case_data\H_Sequence_Limit.csv';

%% Data Loading and Preprocessing
inputConcentrationMatrix = readmatrix(inputDataPath);
pmfFactorScoresMatrix = readmatrix(pmfFactorScoresPath);
pmfFactorLoadingsMatrix = readmatrix(pmfFactorLoadingsPath);

% Extract concentration data (columns 5 to 17 for metals)
concentrationMatrix = inputConcentrationMatrix(:, 5:17);

% Compute PMF-reconstructed concentrations for comparison
pmfReconstructedMatrix = pmfFactorScoresMatrix(:, 2:end) * pmfFactorLoadingsMatrix(:, 2:end)';

% Calculate R² and RMSE for PMF results
pmfR2Scores = zeros(1, numMetals);
pmfRmseValues = zeros(1, numMetals);
for metalIdx = 1:numMetals
    [pmfR2Scores(metalIdx), pmfRmseValues(metalIdx)] = computeR2AndRmse(...
        concentrationMatrix(:, metalIdx), pmfReconstructedMatrix(:, metalIdx));
end

%% Data Normalization (Min-Max Scaling)
normalizedConcentrationMatrix = zeros(numSamples, numMetals);
normalizationParams = zeros(2, numMetals);  % Row 1: min, Row 2: max
for metalIdx = 1:numMetals
    [normalizedCol, colMin, colMax] = minMaxNormalization(concentrationMatrix(:, metalIdx));
    normalizationParams(1, metalIdx) = colMin;
    normalizationParams(2, metalIdx) = colMax;
    normalizedConcentrationMatrix(:, metalIdx) = normalizedCol;
end

%% Load Sequence Constraints for Factor Loadings
sequenceConstraintsTable = readtable(sequenceConstraintsPath);
factorLoadingSequenceLimits = sequenceConstraintsTable{:, 2:end};  % Exclude first column (source names)
factorLoadingSequenceLimits(isnan(factorLoadingSequenceLimits)) = 0;  % Replace NaN with 0

% Define possible factor combinations (e.g., source indices)
factorCombinations = [1, 2, 3, 4; 1, 2, 3, 5];
numCombinations = size(factorCombinations, 1);

%% NMF Execution Across Factor Combinations
resultsCellArray = cell(numCombinations, 10);  % Columns: sourceCombo, factorScores, factorLoadings, errorNorm, constraintSuccessFlag, constraintError, validationFlag, denormErrorNorm, r2Scores, rmseValues
for comboIdx = 1:numCombinations
    fprintf('Processing combination %d of %d\n', comboIdx, numCombinations);
    selectedSources = factorCombinations(comboIdx, :);
    
    % Run multiple replicates to find best factorization
    replicateResults = cell(numReplicates, 5);  % Per replicate: scores, loadings, error, flag, error
    for repIdx = 1:numReplicates
        % Initialize with unconstrained NMF (multiplicative update)
        initialIterativeRecord = cell(maxIterations, 5);
        options = statset('MaxIter', maxIterations, 'Display', 'final');
        [initialFactorScores, initialFactorLoadings, ~, initialIterativeRecord] = nnmfWithConstraints(...
            normalizedConcentrationMatrix, numFactors, [], initialIterativeRecord, ...
            'Replicates', 20, 'Options', options, 'Algorithm', 'mult');
        
        % Apply constraints and refine with ALS
        constrainedSequenceLimits = factorLoadingSequenceLimits(selectedSources, :);
        iterativeRecord = cell(maxIterations, 5);
        options = statset('MaxIter', maxIterations, 'Display', 'final');
        [refinedFactorScores, refinedFactorLoadings, reconstructionError, iterativeRecord] = nnmfWithConstraints(...
            normalizedConcentrationMatrix, numFactors, constrainedSequenceLimits, iterativeRecord, ...
            'W0', initialFactorScores, 'H0', initialFactorLoadings, 'Options', options, 'Algorithm', 'als');
        
        replicateResults{repIdx, :} = findBestFactorization(iterativeRecord, maxIterations);
    end
    
    % Select best replicate
    bestReplicate = findBestFactorization(replicateResults, numReplicates);
    [constraintValidationFlag, constraintValidationError] = validateFactorLoadings(bestReplicate{1, 2}, constrainedSequenceLimits);
    
    % Denormalize reconstructed matrix
    normalizedReconstructed = bestReplicate{1, 1} * bestReplicate{1, 2};
    denormalizedReconstructed = zeros(numSamples, numMetals);
    for metalIdx = 1:numMetals
        denormalizedReconstructed(:, metalIdx) = inverseMinMaxNormalization(...
            normalizedReconstructed(:, metalIdx), normalizationParams(1, metalIdx), normalizationParams(2, metalIdx));
    end
    denormalizedErrorNorm = norm(concentrationMatrix - denormalizedReconstructed, 'fro') / sqrt(numSamples * numMetals);
    
    % Compute R² and RMSE for NMF vs. PMF
    comparisonR2Scores = zeros(2, numMetals);  % Row 1: PMF, Row 2: NMF
    comparisonRmseValues = zeros(2, numMetals);
    comparisonR2Scores(1, :) = pmfR2Scores;
    comparisonRmseValues(1, :) = pmfRmseValues;
    for metalIdx = 1:numMetals
        [comparisonR2Scores(2, metalIdx), comparisonRmseValues(2, matrixIdx)] = computeR2AndRmse(...
            concentrationMatrix(:, metalIdx), denormalizedReconstructed(:, metalIdx));
    end
    
    % Store results
    resultsCellArray{comboIdx, 1} = sequenceConstraintsTable(selectedSources, 1);  % Source names
    resultsCellArray{comboIdx, 2} = bestReplicate{1, 1};  % Factor scores (G)
    resultsCellArray{comboIdx, 3} = bestReplicate{1, 2};  % Factor loadings (F)
    resultsCellArray{comboIdx, 4} = bestReplicate{1, 3};  % Normalized error
    resultsCellArray{comboIdx, 5} = bestReplicate{1, 4};  % Constraint success flag
    resultsCellArray{comboIdx, 6} = bestReplicate{1, 5};  % Constraint error
    resultsCellArray{comboIdx, 7} = constraintValidationFlag;  % Re-validation flag
    resultsCellArray{comboIdx, 8} = denormalizedErrorNorm;  % Denormalized error
    resultsCellArray{comboIdx, 9} = comparisonR2Scores;  % R² comparison
    resultsCellArray{comboIdx, 10} = comparisonRmseValues;  % RMSE comparison
end

% Save results
save('NMF_Source_Apportionment_Results.mat', 'resultsCellArray');

%% Helper Functions
function [r2Score, rmseValue] = computeR2AndRmse(observedValues, predictedValues)
    % Compute R² (coefficient of determination) and RMSE (root mean square error).
    if length(observedValues) ~= length(predictedValues)
        error('Observed and predicted vectors must have the same length.');
    end
    ssRes = sum((observedValues - predictedValues).^2);  % Residual sum of squares
    ssTot = sum((observedValues - mean(observedValues)).^2);  % Total sum of squares
    r2Score = 1 - (ssRes / ssTot);
    rmseValue = sqrt(mean((observedValues - predictedValues).^2));
end

function bestFactorization = findBestFactorization(iterativeRecords, maxRecords)
    % Find the best factorization from iterative records, prioritizing constraint success and minimal error.
    bestFactorization = [];
    minErrorSuccess = inf;
    minErrorFailure = inf;
    foundSuccess = false;
    for recordIdx = 1:maxRecords
        if ~isempty(iterativeRecords{recordIdx, 1})
            currentError = iterativeRecords{recordIdx, 3};
            if iterativeRecords{recordIdx, 4} == 0  % Success flag
                foundSuccess = true;
                if currentError < minErrorSuccess
                    minErrorSuccess = currentError;
                    bestFactorization = iterativeRecords(recordIdx, :);
                end
            elseif ~foundSuccess
                if currentError < minErrorFailure
                    minErrorFailure = currentError;
                    bestFactorization = iterativeRecords(recordIdx, :);
                end
            end
        end
    end
    if foundSuccess
        fprintf('Best successful factorization found.\n');
    else
        fprintf('No successful factorization; using minimal error.\n');
    end
    disp(bestFactorization);
end

function [normalizedVector, minValue, maxValue] = minMaxNormalization(inputVector)
    % Apply min-max normalization to [0, 1] range.
    minValue = min(inputVector);
    maxValue = max(inputVector);
    normalizedVector = (inputVector - minValue) / (maxValue - minValue);
end

function originalVector = inverseMinMaxNormalization(normalizedVector, minValue, maxValue)
    % Inverse min-max normalization to original scale.
    originalVector = normalizedVector * (maxValue - minValue) + minValue;
end

function [factorScores, factorLoadings, errorNorm, iterativeRecord] = nnmfWithConstraints(inputMatrix, numFactors, sequenceLimits, iterativeRecord, varargin)
    % Constrained NMF using alternating least squares (ALS) or multiplicative updates.
    % Supports initial values, replicates, and options. This is a wrapper for the internal single factorization.
    % Inputs:
    %   inputMatrix - Non-negative data matrix (numSamples x numMetals)
    %   numFactors - Number of factors k
    %   sequenceLimits - Constraint matrix for factor loadings (numFactors x numMetals)
    %   iterativeRecord - Cell array for tracking iterations (maxIterations x 5)
    % Varargin: 'algorithm' ('als' or 'mult'), 'W0' (initial factor scores), 'H0' (initial factor loadings),
    %           'Replicates' (int), 'Options' (statset)
    
    if nargin > 4
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
    
    narginchk(4, Inf);
    [numRows, numCols] = size(inputMatrix);
    if ~isscalar(numFactors) || ~isnumeric(numFactors) || numFactors < 1 || numFactors > min(numRows, numCols) || numFactors ~= round(numFactors)
        error('stats:nnmf:BadK');
    end
    
    % Parse optional arguments
    paramNames = {'algorithm', 'w0', 'h0', 'replicates', 'options'};
    defaultValues = {'als', [], [], 1, []};
    [algorithm, initialFactorScores, initialFactorLoadings, numReplicates, options] = internalParseArgs(paramNames, defaultValues, varargin{:});
    
    % Validate algorithm
    algorithm = internalGetParamVal(algorithm, {'mult', 'als'}, 'ALGORITHM');
    isMultiplicative = strncmp('mult', algorithm, length(algorithm));
    
    % Validate matrices
    validateNonNegativeMatrices(inputMatrix, initialFactorScores, initialFactorLoadings, numFactors);
    
    if ~isscalar(numReplicates) || ~isnumeric(numReplicates) || numReplicates < 1 || numReplicates ~= round(numReplicates)
        error('stats:nnmf:BadReplicates');
    end
    
    % Get options
    defaultOpt = statset('nnmf');
    tolX = statget(options, 'TolX', defaultOpt, 'fast');
    tolFun = statget(options, 'TolFun', defaultOpt, 'fast');
    maxIter = statget(options, 'MaxIter', defaultOpt, 'fast');
    dispOpt = statget(options, 'Display', defaultOpt, 'fast');
    
    [~, dispNum] = internalGetParamVal(dispOpt, {'off', 'notify', 'final', 'iter'}, 'Display');
    dispNum = dispNum - 1;
    
    % Parallel processing setup (simplified; assumes no pool for base case)
    useParallel = false;  % Extend if needed
    usePool = false;
    
    % Special case for full rank
    if isempty(initialFactorScores) && isempty(initialFactorLoadings)
        if numFactors == numCols
            initialFactorScores = inputMatrix;
            initialFactorLoadings = eye(numFactors);
        elseif numFactors == numRows
            initialFactorScores = eye(numFactors);
            initialFactorLoadings = inputMatrix;
        end
    end
    
    % Define loop body for replicates
    loopBodyFcn = @(repIdx, rngStream) singleNmfReplicate(repIdx, rngStream, inputMatrix, numFactors, sequenceLimits, ...
        initialFactorScores, initialFactorLoadings, isMultiplicative, maxIter, tolFun, tolX, dispNum, usePool, iterativeRecord);
    
    % Run replicates and find best
    try
        [bestNorm, ~, bestScores, bestLoadings] = runReplicates(numReplicates, loopBodyFcn, useParallel);
    catch ME
        rethrow(ME);
    end
    
    errorNorm = bestNorm;
    factorScores = bestScores;
    factorLoadings = bestLoadings;
    
    if dispNum > 1
        fprintf('Final RMS Residual: %g\n', bestNorm);
    end
    
    if bestNorm == Inf
        error('stats:nnmf:NoSolution');
    end
    
    % Normalize and order factors (standard NMF post-processing)
    hLen = sqrt(sum(bestLoadings.^2, 2));
    if any(hLen == 0)
        warning('stats:nnmf:LowRank', 'Low rank: %d of %d factors', numFactors - sum(hLen == 0), numFactors);
        hLen(hLen == 0) = 1;
    end
    factorScores = bsxfun(@times, factorScores, hLen');
    factorLoadings = bsxfun(@times, factorLoadings, 1 ./ hLen);
    
    [~, sortIdx] = sort(sum(factorScores.^2, 1), 'descend');
    factorScores = factorScores(:, sortIdx);
    factorLoadings = factorLoadings(sortIdx, :);
    
    function [repNorm, repIdx, repScores, repLoadings] = singleNmfReplicate(repIdx, rngStream, inputMatrix, numFactors, sequenceLimits, ...
            initialScores, initialLoadings, isMult, maxIter, tolFun, tolX, dispNum, usePool, iterativeRecord)
        % Single replicate NMF
        if isempty(rngStream)
            rngStream = RandStream.getGlobalStream;
        end
        
        % Initialize
        currentScores = initialScores;
        if isempty(currentScores)
            currentScores = rand(rngStream, numRows, numFactors);
        end
        currentLoadings = initialLoadings;
        if isempty(currentLoadings)
            currentLoadings = rand(rngStream, numFactors, numCols);
        end
        
        % Perform factorization
        [currentScores, currentLoadings, currentNorm, iterativeRecord] = singleNmfFactorization(...
            inputMatrix, currentScores, currentLoadings, isMult, maxIter, tolFun, tolX, dispNum, repIdx, usePool, sequenceLimits, iterativeRecord);
        
        repNorm = currentNorm;
        repIdx = repIdx;
        repScores = currentScores;
        repLoadings = currentLoadings;
    end
    
    function [bestNorm, bestRepIdx, bestScores, bestLoadings] = runReplicates(numReps, loopBody, usePar)
        % Run and reduce replicates (simplified serial version)
        bestNorm = inf;
        bestRepIdx = 1;
        bestScores = [];
        bestLoadings = [];
        for rep = 1:numReps
            [repNorm, repIdx, repScores, repLoadings] = loopBody(rep, []);
            if repNorm < bestNorm
                bestNorm = repNorm;
                bestRepIdx = repIdx;
                bestScores = repScores;
                bestLoadings = repLoadings;
            end
        end
    end
end

function [scores, loadings, normError, iterativeRecord] = singleNmfFactorization(inputMatrix, initialScores, initialLoadings, ...
        isMultiplicative, maxIter, tolFun, tolX, dispNum, repNum, usePool, sequenceLimits, iterativeRecord)
    % Core single NMF iteration with constraint clamping on loadings.
    numElements = numel(inputMatrix);
    sqrtEps = sqrt(eps);
    
    % Display setup
    if dispNum > 1
        if usePool
            labIdx = 0;  % Placeholder
            dispFmt = '%8d\t%8d\t%8d\t%14g\t%14g\n';
        else
            dispFmt = '%7d\t%8d\t%12g\t%12g\n';
        end
    end
    
    prevNormError = 0;
    currentScores = initialScores;
    currentLoadings = initialLoadings;
    
    for iter = 1:maxIter
        if isMultiplicative
            % Multiplicative update
            numerator = currentScores' * inputMatrix;
            currentLoadings = max(0, currentLoadings .* (numerator ./ ((currentScores' * currentScores) * currentLoadings + eps(numerator))));
            [currentLoadings, successFlag, constraintError] = clampFactorLoadings(currentLoadings, sequenceLimits);
            numerator = inputMatrix * currentLoadings';
            currentScores = max(0, currentScores .* (numerator ./ (currentScores * (currentLoadings * currentLoadings') + eps(numerator))));
        else
            % ALS update
            currentLoadings = max(0, currentScores \ inputMatrix);
            [currentLoadings, successFlag, constraintError] = clampFactorLoadings(currentLoadings, sequenceLimits);
            currentScores = max(0, inputMatrix / currentLoadings);
        end
        
        % Compute error and delta
        residual = inputMatrix - currentScores * currentLoadings;
        normError = sqrt(sum(sum(residual.^2)) / numElements);
        deltaScores = max(max(abs(currentScores - initialScores) / (sqrtEps + max(max(abs(initialScores))))));
        deltaLoadings = max(max(abs(currentLoadings - initialLoadings) / (sqrtEps + max(max(abs(initialLoadings))))));
        delta = max(deltaScores, deltaLoadings);
        
        % Record iteration
        iterativeRecord{iter, 1} = currentScores;
        iterativeRecord{iter, 2} = currentLoadings;
        iterativeRecord{iter, 3} = normError;
        iterativeRecord{iter, 4} = successFlag;
        iterativeRecord{iter, 5} = constraintError;
        
        % Convergence check
        if iter > 1
            if delta <= tolX
                break;
            elseif prevNormError - normError <= tolFun * max(1, prevNormError)
                break;
            elseif iter == maxIter
                break;
            end
        end
        
        if dispNum > 2  % 'iter'
            if usePool
                fprintf(dispFmt, labIdx, repNum, iter, normError, delta);
            else
                fprintf(dispFmt, repNum, iter, normError, delta);
            end
        end
        
        % Update for next iteration
        prevNormError = normError;
        initialScores = currentScores;
        initialLoadings = currentLoadings;
    end
    
    if dispNum > 1
        if usePool
            fprintf(dispFmt, labIdx, repNum, iter, normError, delta);
        else
            fprintf(dispFmt, repNum, iter, normError, delta);
        end
    end
    
    scores = currentScores;
    loadings = currentLoadings;
end

function [clampedLoadings, successFlag, errorMatrix] = clampFactorLoadings(factorLoadings, sequenceLimits)
    % Clamp factor loadings to satisfy sequence constraints using constrained optimization.
    numFactorsLocal = size(factorLoadings, 1);
    clampedLoadings = factorLoadings;
    
    if isempty(sequenceLimits)
        successFlag = NaN;
        errorMatrix = NaN;
        return;
    end
    
    % Identify columns with non-trivial constraints
    nonZeroCols = find(any(sequenceLimits ~= 0, 1));
    numNonZeroCols = length(nonZeroCols);
    
    successFlag = 0;
    errorMatrix = zeros(size(factorLoadings));
    
    for colIdx = 1:numNonZeroCols
        colNum = nonZeroCols(colIdx);
        if all(sequenceLimits(:, colNum) == 1)
            continue;  % No constraint
        end
        
        colSum = sum(factorLoadings(:, colNum));
        objectiveFcn = @(x) sum(abs(x - factorLoadings(:, colNum)));  % Minimize deviation
        
        % Generate inequality constraints A*x <= b (from sequence)
        [constraintMatrixA, numConstraints] = generateSequenceConstraints(sequenceLimits(:, colNum), numFactorsLocal);
        bVec = zeros(numConstraints, 1);
        aEq = ones(1, numFactorsLocal);  % Sum equality
        bEq = colSum;
        lb = zeros(numFactorsLocal, 1);
        ub = colSum * ones(numFactorsLocal, 1);
        
        % Optimize with fmincon (try up to 10 initial points)
        converged = false;
        for initTry = 1:10
            if initTry == 1
                x0 = factorLoadings(:, colNum);
            else
                x0 = unifrnd(0, colSum, numFactorsLocal, 1);
            end
            [xOpt, ~, exitFlag] = fmincon(objectiveFcn, x0, constraintMatrixA, bVec, aEq, bEq, lb, ub);
            if exitFlag > 0
                clampedLoadings(:, colNum) = xOpt;
                converged = true;
                break;
            end
        end
        if ~converged
            warning('Constraint optimization failed for column %d', colNum);
        end
    end
    
    % Re-validate
    [successFlag, errorMatrix] = validateFactorLoadings(clampedLoadings, sequenceLimits);
    if successFlag ~= 0
        warning('Post-clamping validation failed');
    end
end

function [validationFlag, errorMatrix] = validateFactorLoadings(factorLoadings, sequenceLimits)
    % Validate if factor loadings satisfy sequence constraints.
    numFactorsLocal = size(factorLoadings, 1);
    validationFlag = 0;
    errorMatrix = zeros(size(factorLoadings));
    
    if isempty(sequenceLimits)
        validationFlag = NaN;
        errorMatrix = NaN;
        return;
    end
    
    nonZeroCols = find(any(sequenceLimits ~= 0, 1));
    numNonZeroCols = length(nonZeroCols);
    
    for colIdx = 1:numNonZeroCols
        colNum = nonZeroCols(colIdx);
        if all(sequenceLimits(:, colNum) == 1)
            continue;
        end
        [colFlag, colError] = validateSequenceConstraints(factorLoadings(:, colNum), sequenceLimits(:, colNum), numFactorsLocal);
        validationFlag = validationFlag + colFlag;
        errorMatrix(:, colNum) = colError;
    end
end

function [colFlag, colError] = validateSequenceConstraints(colLoadings, sequenceLimits, numFactorsLocal)
    % Validate a single column against sequence limits.
    [isPositions, outPositions, maxNum] = findNonSpecificPositions(sequenceLimits, numFactorsLocal);
    colFlag = 0;
    colError = zeros(numFactorsLocal, 1);
    
    for groupIdx = 1:maxNum
        numIs = length(isPositions{groupIdx});
        numOut = length(outPositions{groupIdx});
        for isIdx = 1:numIs
            for outIdx = 1:numOut
                if colLoadings(isPositions{groupIdx}(isIdx)) < colLoadings(outPositions{groupIdx}(outIdx))
                    colFlag = colFlag + 1;
                    colError(isPositions{groupIdx}(isIdx)) = colError(isPositions{groupIdx}(isIdx)) + colLoadings(isPositions{groupIdx}(isIdx));
                    colError(outPositions{groupIdx}(outIdx)) = colError(outPositions{groupIdx}(outIdx)) + colLoadings(outPositions{groupIdx}(outIdx));
                end
            end
        end
    end
end

function [constraintA, numCons] = generateSequenceConstraints(sequenceLimits, numFactorsLocal)
    % Generate inequality constraints A from sequence limits.
    [isPositions, outPositions, maxNum] = findNonSpecificPositions(sequenceLimits, numFactorsLocal);
    constraintA = zeros(numFactorsLocal * numFactorsLocal, numFactorsLocal);
    numCons = 0;
    
    for groupIdx = 1:maxNum
        numIs = length(isPositions{groupIdx});
        numOut = length(outPositions{groupIdx});
        for isIdx = 1:numIs
            for outIdx = 1:numOut
                numCons = numCons + 1;
                constraintA(numCons, isPositions{groupIdx}(isIdx)) = -1;
                constraintA(numCons, outPositions{groupIdx}(outIdx)) = 1;
            end
        end
    end
    constraintA = constraintA(1:numCons, :);
end

function [isPositions, outPositions, maxNum] = findNonSpecificPositions(sequenceLimits, numFactorsLocal)
    % Identify positions for "is less than" constraints based on sequence values.
    isPositions = cell(1, numFactorsLocal);
    outPositions = cell(1, numFactorsLocal);
    mask = true(size(sequenceLimits));
    maxNum = 0;
    
    for val = 1:numFactorsLocal
        mask = mask & (sequenceLimits ~= val);
        isPos = find(sequenceLimits == val);
        if isempty(isPos)
            maxNum = val - 1;
            break;
        end
        outPos = find(mask);
        isPositions{val} = isPos;
        outPositions{val} = outPos;
        maxNum = val;
    end
end

function validateNonNegativeMatrices(inputMat, scoresMat, loadingsMat, numFactorsLocal)
    % Validate input matrices are non-negative and properly sized.
    if ~ismatrix(inputMat) || ~isnumeric(inputMat) || ~isreal(inputMat) || any(any(~isfinite(inputMat)))
        error('stats:nnmf:BadA');
    end
    [numR, numC] = size(inputMat);
    
    if ~isempty(scoresMat)
        if ~ismatrix(scoresMat) || ~isnumeric(scoresMat) || ~isreal(scoresMat) || any(any(scoresMat < 0)) || any(any(~isfinite(scoresMat)))
            error('stats:nnmf:BadWNegativeValues');
        elseif ~isequal(size(scoresMat), [numR, numFactorsLocal])
            error('stats:nnmf:BadWSizeIsWrong', num2str(numR), num2str(numFactorsLocal));
        end
    end
    
    if ~isempty(loadingsMat)
        if ~ismatrix(loadingsMat) || ~isnumeric(loadingsMat) || ~isreal(loadingsMat) || any(any(loadingsMat < 0)) || any(any(~isfinite(loadingsMat)))
            error('stats:nnmf:BadHNegativeValues');
        elseif ~isequal(size(loadingsMat), [numFactorsLocal, numC])
            error('stats:nnmf:BadHSizeIsWrong', num2str(numFactorsLocal), num2str(numC));
        end
    end
end

function [params, values] = internalParseArgs(paramNames, defaultValues, varargin)
    % Internal argument parser (simplified).
    params = paramNames;
    values = defaultValues;
    % Parse varargin into params/values (MATLAB internal style; implement as needed).
end

function paramVal = internalGetParamVal(val, validVals, paramName)
    % Internal parameter validation.
    for i = 1:length(validVals)
        if strcmpi(val, validVals{i})
            paramVal = validVals{i};
            return;
        end
    end
    error('stats:nnmf:Bad%s', paramName);
end