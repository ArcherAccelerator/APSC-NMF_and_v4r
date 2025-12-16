clear; clc;
%% NMF Uncertainty Analysis for Source Apportionment
% This script quantifies uncertainty in NMF factor scores and loadings 
% via Monte Carlo replicates, computing mean and standard deviation.
% Builds on NMF_Source_Apportionment.m with repeated runs for uncertainty estimation.
% Author: Yuqinqin (Modified for clarity and professionalism, 2025)
% Date: November 2, 2024 (Original); December 16, 2025 (Revised)

%% Configuration Parameters
numSamples = 95;
numMetals = 13;
numFactors = 4;
maxIterations = 1000;
numReplicatesPerRun = 200;
numUncertaintyRuns = 10;  % Monte Carlo runs for uncertainty

% File paths
inputDataPath = '..\Case_data\SourceApportionmentData.xls';
sequenceConstraintsPath = 'H_Sequence_Limit.csv';

%% Data Loading and Normalization
inputConcentrationMatrix = readmatrix(inputDataPath);
concentrationMatrix = inputConcentrationMatrix(:, 5:17);  % Extract metals
normalizedConcentrationMatrix = zeros(numSamples, numMetals);
normalizationParams = zeros(2, numMetals);
for metalIdx = 1:numMetals
    [normCol, colMin, colMax] = minMaxNormalization(concentrationMatrix(:, metalIdx));
    normalizationParams(1, metalIdx) = colMin;
    normalizationParams(2, metalIdx) = colMax;
    normalizedConcentrationMatrix(:, metalIdx) = normCol;
end

%% Load Constraints
sequenceConstraintsTable = readtable(sequenceConstraintsPath);
factorLoadingSequenceLimits = sequenceConstraintsTable{:, 2:end};
factorLoadingSequenceLimits(isnan(factorLoadingSequenceLimits)) = 0;

% Selected factor combination (e.g., sources 1,2,4,5 for this file)
selectedSources = [1, 2, 4, 5];

%% Monte Carlo Uncertainty Analysis
uncertaintyResultsCell = cell(numUncertaintyRuns, 8);  % Per run: sourceCombo, scores, loadings, normError, successFlag, error, validationFlag, denormError
factorScoresCollection = cell(numUncertaintyRuns, 1);
factorLoadingsCollection = cell(numUncertaintyRuns, 1);
denormalizedErrors = zeros(numUncertaintyRuns, 1);

for runIdx = 1:numUncertaintyRuns
    fprintf('Processing uncertainty run %d of %d\n', runIdx, numUncertaintyRuns);
    constrainedSequenceLimits = factorLoadingSequenceLimits(selectedSources, :);
    
    % Run replicates to find best per run
    replicateResults = cell(numReplicatesPerRun, 5);
    for repIdx = 1:numReplicatesPerRun
        % Unconstrained initialization (multiplicative)
        initialIterativeRecord = cell(maxIterations, 5);
        options = statset('MaxIter', maxIterations, 'Display', 'final');
        [initialFactorScores, initialFactorLoadings, ~, initialIterativeRecord] = nnmfWithConstraints(...
            normalizedConcentrationMatrix, numFactors, [], initialIterativeRecord, ...
            'Replicates', 20, 'Options', options, 'Algorithm', 'mult');
        
        % Constrained refinement (ALS)
        iterativeRecord = cell(maxIterations, 5);
        options = statset('MaxIter', maxIterations, 'Display', 'final');
        [refinedFactorScores, refinedFactorLoadings, reconstructionError, iterativeRecord] = nnmfWithConstraints(...
            normalizedConcentrationMatrix, numFactors, constrainedSequenceLimits, iterativeRecord, ...
            'W0', initialFactorScores, 'H0', initialFactorLoadings, 'Options', options, 'Algorithm', 'als');
        
        replicateResults{repIdx, :} = findBestFactorization(iterativeRecord, maxIterations);
    end
    
    % Best per run
    bestReplicate = findBestFactorization(replicateResults, numReplicatesPerRun);
    [constraintValidationFlag, constraintValidationError] = validateFactorLoadings(bestReplicate{1, 2}, constrainedSequenceLimits);
    
    % Denormalize for error computation
    normalizedReconstructed = bestReplicate{1, 1} * bestReplicate{1, 2};
    denormalizedReconstructed = zeros(numSamples, numMetals);
    for metalIdx = 1:numMetals
        denormalizedReconstructed(:, metalIdx) = inverseMinMaxNormalization(...
            normalizedReconstructed(:, metalIdx), normalizationParams(1, metalIdx), normalizationParams(2, metalIdx));
    end
    currentDenormError = norm(concentrationMatrix - denormalizedReconstructed, 'fro') / sqrt(numSamples * numMetals);
    
    % Store
    uncertaintyResultsCell{runIdx, 1} = sequenceConstraintsTable(selectedSources, 1);  % Source names
    uncertaintyResultsCell{runIdx, 2} = bestReplicate{1, 1};  % Factor scores
    uncertaintyResultsCell{runIdx, 3} = bestReplicate{1, 2};  % Factor loadings
    uncertaintyResultsCell{runIdx, 4} = bestReplicate{1, 3};  % Normalized error
    uncertaintyResultsCell{runIdx, 5} = bestReplicate{1, 4};  % Success flag
    uncertaintyResultsCell{runIdx, 6} = bestReplicate{1, 5};  % Constraint error
    uncertaintyResultsCell{runIdx, 7} = constraintValidationFlag;  % Validation flag
    uncertaintyResultsCell{runIdx, 8} = currentDenormError;  % Denormalized error
    
    factorScoresCollection{runIdx} = bestReplicate{1, 1};
    factorLoadingsCollection{runIdx} = bestReplicate{1, 2};
    denormalizedErrors(runIdx) = currentDenormError;
end

%% Compute Uncertainty Metrics
[meanFactorScores, stdFactorScores] = computeUncertaintyMetrics(factorScoresCollection);
[meanFactorLoadings, stdFactorLoadings] = computeUncertaintyMetrics(factorLoadingsCollection);

% Save results
save('NMF_Uncertainty_Metrics_5.mat', 'meanFactorScores', 'meanFactorLoadings', 'stdFactorScores', 'stdFactorLoadings', 'denormalizedErrors', 'uncertaintyResultsCell');

%% Helper Functions
function [meanMatrix, stdMatrix] = computeUncertaintyMetrics(matrixCollection)
    % Compute element-wise mean and standard deviation across Monte Carlo runs.
    numRuns = length(matrixCollection);
    [numRows, numCols] = size(matrixCollection{1});
    stackedMatrices = zeros(numRows, numCols, numRuns);
    for runIdx = 1:numRuns
        stackedMatrices(:, :, runIdx) = matrixCollection{runIdx};
    end
    meanMatrix = mean(stackedMatrices, 3);
    stdMatrix = std(stackedMatrices, 0, 3);
end

% Reuse all helper functions from NMF_Source_Apportionment.m:
% - computeR2AndRmse (not used here, but available)
% - findBestFactorization
% - minMaxNormalization
% - inverseMinMaxNormalization
% - nnmfWithConstraints (full implementation)
% - singleNmfFactorization
% - clampFactorLoadings
% - validateFactorLoadings
% - validateSequenceConstraints
% - generateSequenceConstraints
% - findNonSpecificPositions
% - validateNonNegativeMatrices
% - internalParseArgs
% - internalGetParamVal