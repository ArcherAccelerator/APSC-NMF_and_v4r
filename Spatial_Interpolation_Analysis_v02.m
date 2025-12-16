clear; clc; close all;
%% Spatial Interpolation Analysis for Heavy Metal Concentrations
% This script performs comparative spatial interpolation (OK, IDW, v4, v4r) 
% on heavy metal data with cross-validation using random sampling.
% Outputs: Visualization figures for interpolation performance.
% Author: [Your Name] (Modified for clarity and professionalism, 2025)
% Date: April 11, 2025 (Original); December 16, 2025 (Revised)

%% Configuration Parameters
dataPath = '..\Case_data\SamplingPoint_HeavyMetalsConcentrations.xlsx';
classifiedDataPath = '..\Case_data\SamplingData_InteriorPoints.xlsx';
shapefilePath = '..\Case_data\Q4.shp';  % Study area boundary
outputDir = '..\Figure';
addpath('..\ooDACE-1.4\ooDACE');  % Kriging toolbox

% Interpolation methods and metals
interpolationMethods = {'OK', 'IDW', 'v4', 'v4r'};
heavyMetals = {'Cd', 'Hg', 'As', 'Pb', 'Cr', 'Cu', 'Ni', 'Zn', 'Tl', 'Be', 'Se', 'Sb'};
sampleSizes = [5, 10, 15, 20];  % Cross-validation sample sizes
gridResolution = 200;  % X-axis sampling points

%% Data Loading
concentrationData = readmatrix(dataPath);
sampleLongitudes = concentrationData(:, 1);
sampleLatitudes = concentrationData(:, 2);
trueConcentrationsMatrix = concentrationData(:, 3:end);  % Columns: metals

% Load point classifications (e.g., interior vs. boundary)
classifiedData = readmatrix(classifiedDataPath);
pointClassifications = classifiedData(:, 1);

%% Load Study Area Boundary (Shapefile)
boundaryShapes = shaperead(shapefilePath);
boundaryLongitudes = []; boundaryLatitudes = [];
for shapeIdx = 1:length(boundaryShapes)
    boundaryLongitudes = [boundaryLongitudes; boundaryShapes(shapeIdx).X(:); NaN];
    boundaryLatitudes = [boundaryLatitudes; boundaryShapes(shapeIdx).Y(:); NaN];
end
validBoundaryIdx = ~(isnan(boundaryLongitudes) | isnan(boundaryLatitudes));
boundaryLongitudes = boundaryLongitudes(validBoundaryIdx);
boundaryLatitudes = boundaryLatitudes(validBoundaryIdx);

%% Define Plotting Extent (with Padding)
minLon = min(boundaryLongitudes); maxLon = max(boundaryLongitudes);
minLat = min(boundaryLatitudes); maxLat = max(boundaryLatitudes);
centerLon = (minLon + maxLon) / 2;
centerLat = (minLat + maxLat) / 2;
maxExtent = max(maxLon - minLon, maxLat - minLat);
paddingRatio = 0.1;
halfPaddedExtent = (maxExtent / 2) * (1 + paddingRatio);
plotMinLon = centerLon - halfPaddedExtent;
plotMaxLon = centerLon + halfPaddedExtent;
plotMinLat = centerLat - halfPaddedExtent;
plotMaxLat = centerLat + halfPaddedExtent;

%% Generate Uniform Interpolation Grid
[lonGridRange, latGridRange, gridDensity] = generateEquispacedGrid(...
    concentrationData, gridResolution, minLon, maxLon, minLat, maxLat);
[queryLongitudes, queryLatitudes] = meshgrid(lonGridRange, latGridRange);

%% Create Mask for Study Area Clipping
studyAreaMask = false(size(queryLongitudes));
inBoundaryMask = false(size(sampleLatitudes));
for shapeIdx = 1:length(boundaryShapes)
    polyLons = boundaryShapes(shapeIdx).X;
    polyLats = boundaryShapes(shapeIdx).Y;
    inPolyGrid = inpolygon(queryLongitudes, queryLatitudes, polyLons, polyLats);
    studyAreaMask = studyAreaMask | inPolyGrid;
    inBoundaryMask = inBoundaryMask | inpolygon(sampleLongitudes, sampleLatitudes, polyLons, polyLats);
end

%% Main Cross-Validation Loop
numSampleSizes = length(sampleSizes);
numMetals = length(heavyMetals);
for sizeIdx = 1:numSampleSizes
    targetSampleSize = sampleSizes(sizeIdx);
    for metalIdx = 1:numMetals
        trueConcentrations = trueConcentrationsMatrix(:, metalIdx);
        validDataIdx = ~isnan(trueConcentrations);
        validLongitudes = sampleLongitudes(validDataIdx);
        validLatitudes = sampleLatitudes(validDataIdx);
        validConcentrations = trueConcentrations(validDataIdx);
        validClassifications = pointClassifications(validDataIdx);
        
        numValidPoints = length(validConcentrations);
        [interiorPointIndices] = identifyInteriorPoints([validClassifications, validLongitudes, validLatitudes]);
        [samplingResults] = generateRandomSamples(targetSampleSize, interiorPointIndices);
        
        for sampleRunIdx = 1:length(samplingResults)
            fprintf('Processing: Sample size %d, Metal %s, Run %d/%d\n', ...
                targetSampleSize, heavyMetals{metalIdx}, sampleRunIdx, length(samplingResults));
            
            testIndices = samplingResults{sampleRunIdx};
            trainingIndices = setdiff(1:numValidPoints, testIndices);
            
            trainingLongitudes = validLongitudes(trainingIndices);
            trainingLatitudes = validLatitudes(trainingIndices);
            trainingConcentrations = validConcentrations(trainingIndices);
            testLongitudes = validLongitudes(testIndices);
            testLatitudes = validLatitudes(testIndices);
            testConcentrations = validConcentrations(testIndices);
            
            % Perform interpolations
            ordinaryKrigingGrid = ordinaryKriging(trainingLongitudes, trainingLatitudes, trainingConcentrations, queryLongitudes, queryLatitudes);
            inverseDistanceWeightedGrid = inverseDistanceWeighting(trainingLongitudes, trainingLatitudes, trainingConcentrations, queryLongitudes, queryLatitudes);
            linearV4Grid = scatteredInterpolantGrid(trainingLongitudes, trainingLatitudes, trainingConcentrations, queryLongitudes, queryLatitudes, 'linear', 'v4');
            refinedV4RGrid = iterativeNegativeRefinement(trainingLongitudes, trainingLatitudes, trainingConcentrations, queryLongitudes, queryLatitudes);
            
            % Visualize results
            interpolatedGrids = {ordinaryKrigingGrid, inverseDistanceWeightedGrid, linearV4Grid, refinedV4RGrid};
            methodNames = {'OK', 'IDW', 'v4', 'v4r'};
            figure('Position', [100, 100, 2000, 1600], 'Visible', 'off');
            for methodIdx = 1:length(methodNames)
                subplot(2, 2, methodIdx); hold on;
                interpolatedGrid = interpolatedGrids{methodIdx};
                interpolatedGrid(~studyAreaMask) = NaN;
                [trimmedGrid, trimmedQueryLons, trimmedQueryLats] = trimNaNBoundaries(interpolatedGrid, queryLongitudes, queryLatitudes);
                [trimmedGrid, trimmedQueryLons, trimmedQueryLats] = enforceSquareAspect(trimmedGrid, trimmedQueryLons, trimmedQueryLats, gridDensity);
                
                minVal = min(trimmedGrid(:), [], 'omitnan');
                maxVal = max(trimmedGrid(:), [], 'omitnan');
                h = pcolor(trimmedQueryLons, trimmedQueryLats, trimmedGrid);
                set(h, 'EdgeColor', 'none'); shading interp; colormap(jet);
                c = colorbar; c.Label.FontSize = 18;
                caxis([min(validConcentrations) max(validConcentrations)]);
                
                % Plot training/test points
                inTrainingMask = inpolygon(trainingLongitudes, trainingLatitudes, boundaryLongitudes, boundaryLatitudes);
                scatter(trainingLongitudes(inTrainingMask), trainingLatitudes(inTrainingMask), 80, trainingConcentrations(inTrainingMask), 'filled', 'o', 'MarkerEdgeColor', 'k');
                scatter(testLongitudes, testLatitudes, 80, testConcentrations, 'filled', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
                
                title([methodNames{methodIdx}, sprintf(' (min=%.4f, max=%.4f)']], 'FontSize', 22);
                xlabel('Longitude', 'FontSize', 18); ylabel('Latitude', 'FontSize', 18);
                ax = gca; ax.FontSize = 16; axis equal tight;
                xlim([plotMinLon plotMaxLon]); ylim([plotMinLat plotMaxLat]);
            end
            sgtitle([heavyMetals{metalIdx}, sprintf(' (n=%d, sample %d/%d)')], 'FontSize', 24);
            
            % Save figure
            if ~exist(outputDir, 'dir'), mkdir(outputDir); end
            saveas(gcf, fullfile(outputDir, [heavyMetals{metalIdx}, '_n', num2str(targetSampleSize), '_sample', num2str(sampleRunIdx), '.png']));
            close(gcf);
        end
    end
end

%% Helper Functions
function Z = ordinaryKriging(x, y, z, xq, yq)
    % Ordinary Kriging interpolation using DACE toolbox.
    xyCoords = [x, y];
    thetaInit = [10, 10]; lob = [0.1, 0.1]; upb = [20, 20];
    dmodel = dacefit(xyCoords, z, 'regpoly0', 'corrgauss', thetaInit, lob, upb);
    Z = zeros(size(xq));
    for i = 1:numel(xq)
        [Z(i), ~] = predictor([xq(i), yq(i)], dmodel);
    end
end

function Z = inverseDistanceWeighting(x, y, z, xq, yq)
    % Inverse Distance Weighting (IDW) interpolation.
    Z = zeros(size(xq)); power = 2;
    for i = 1:numel(xq)
        distances = sqrt((x - xq(i)).^2 + (y - yq(i)).^2);
        distances(distances == 0) = eps;
        weights = 1 ./ (distances .^ power);
        Z(i) = sum(weights .* z) / sum(weights);
    end
end

function Z = scatteredInterpolantGrid(x, y, z, xq, yq, method, extrapMethod)
    % MATLAB's scatteredInterpolant for 'v4' (linear) or similar.
    F = scatteredInterpolant(x, y, z, method, extrapMethod);
    Z = F(xq, yq);
end

function Z = iterativeNegativeRefinement(x, y, z, xq, yq, maxIterations)
    % Iterative v4 with negative value refinement (add points to eliminate negatives).
    if nargin < 6, maxIterations = 20; end
    currentX = x; currentY = y; currentZ = z;
    Z = scatteredInterpolantGrid(currentX, currentY, currentZ, xq, yq, 'linear', 'v4');
    for iter = 1:maxIterations
        if min(Z(:)) < 0
            [minVal, linIdx] = min(Z(:));
            [row, col] = ind2sub(size(Z), linIdx);
            currentX = [currentX; xq(row, col)];
            currentY = [currentY; yq(row, col)];
            currentZ = [currentZ; min(currentZ)];  % Clamp to min observed
            Z = scatteredInterpolantGrid(currentX, currentY, currentZ, xq, yq, 'linear', 'v4');
        else
            break;
        end
    end
end

function [samples] = generateRandomSamples(sampleSize, availableIndices)
    % Generate random subsets of indices for cross-validation.
    numAvailable = length(availableIndices);
    numFullSets = floor(numAvailable / sampleSize);
    samples = cell(numFullSets, 1);
    permutedIndices = randperm(numAvailable);
    for setIdx = 1:numFullSets
        startIdx = (setIdx - 1) * sampleSize + 1;
        endIdx = setIdx * sampleSize;
        samples{setIdx} = availableIndices(permutedIndices(startIdx:endIdx));
    end
    remaining = mod(numAvailable, sampleSize);
    if remaining > 0
        remainingIdx = permutedIndices(numFullSets * sampleSize + 1 : end);
        samples{end + 1} = availableIndices(remainingIdx);
    end
end

function [interiorIndices] = identifyInteriorPoints(pointData)
    % Identify points not on convex hull (interior points) by classification.
    uniqueClasses = unique(pointData(:, 1));
    interiorIndices = [];
    for classIdx = 4:length(uniqueClasses)  % Focus on class 4 (interior)
        classData = pointData(pointData(:, 1) == uniqueClasses(classIdx), :);
        pointCoords = classData(:, 2:3);
        hullIndices = convhull(pointCoords(:, 1), pointCoords(:, 2));
        hullPoints = unique(hullIndices);
        classRows = find(pointData(:, 1) == uniqueClasses(classIdx));
        notOnHull = setdiff(1:size(pointCoords, 1), hullPoints);
        interiorIndices = [interiorIndices; classRows(notOnHull)];
    end
end

function [lonRange, latRange, density] = generateEquispacedGrid(data, xResolution, lonMin, lonMax, latMin, latMax)
    % Generate equispaced grid based on data extent and boundary.
    xLimits = [lonMin, lonMax]; yLimits = [latMin, latMax];
    finalMinLon = min(xLimits(1), min(data(:, 1)));
    finalMaxLon = max(xLimits(2), max(data(:, 1)));
    finalMinLat = min(yLimits(1), min(data(:, 2)));
    finalMaxLat = max(yLimits(2), max(data(:, 2)));
    lonLength = finalMaxLon - finalMinLon;
    latLength = finalMaxLat - finalMinLat;
    density = lonLength / xResolution;
    yResolution = ceil(latLength / density);
    latLengthAdjusted = density * yResolution;
    finalMaxLatAdjusted = finalMinLat + latLengthAdjusted;
    lonRange = linspace(finalMinLon, finalMaxLon, xResolution);
    latRange = linspace(finalMinLat, finalMaxLatAdjusted, yResolution);
end

function [trimmedMat, trimmedX, trimmedY] = trimNaNBoundaries(mat, x, y)
    % Trim NaN edges from interpolated matrix.
    while all(isnan(mat(1, :))), mat(1, :) = []; x(1, :) = []; y(1, :) = []; end
    while all(isnan(mat(end, :))), mat(end, :) = []; x(end, :) = []; y(end, :) = []; end
    while all(isnan(mat(:, 1))), mat(:, 1) = []; x(:, 1) = []; y(:, 1) = []; end
    while all(isnan(mat(:, end))), mat(:, end) = []; x(:, end) = []; y(:, end) = []; end
    trimmedMat = mat; trimmedX = x; trimmedY = y;
end

function [Z, xq, yq] = enforceSquareAspect(Z, xq, yq, density)
    % Pad matrix to enforce square aspect ratio for plotting.
    [rows, cols] = size(Z);
    if rows > cols
        diff = rows - cols; q1 = floor(diff / 2); q2 = diff - q1;
        [Z, xq, yq] = padWithNaNColumns(Z, xq, yq, q1, q2, density);
    else
        diff = cols - rows; q1 = floor(diff / 2); q2 = diff - q1;
        [Z, xq, yq] = padWithNaNRows(Z, xq, yq, q1, q2, density);
    end
end

function [Z, xq, yq] = padWithNaNColumns(Z, xq, yq, leftCols, rightCols, density)
    [m, ~] = size(Z);
    leftPad = NaN(m, leftCols); rightPad = NaN(m, rightCols);
    Z = [leftPad, Z, rightPad];
    yq = [yq(:, 1:leftCols), yq, yq(:, 1:rightCols)];
    leftXPad = ones(m, leftCols) * (xq(1, 1) - density * (leftCols:-1:1)');
    rightXPad = ones(m, rightCols) * (xq(end, end) + density * (1:rightCols));
    xq = [leftXPad, xq, rightXPad];
end

function [Z, xq, yq] = padWithNaNRows(Z, xq, yq, topRows, bottomRows, density)
    [~, n] = size(Z);
    topPad = NaN(topRows, n); bottomPad = NaN(bottomRows, n);
    Z = [topPad; Z; bottomPad];
    xq = [xq(1:topRows, :); xq; xq(1:bottomRows, :)];
    topYPad = (yq(1, 1) - density * (topRows:-1:1)') * ones(topRows, n);
    bottomYPad = (yq(end, end) + density * (1:bottomRows)) * ones(bottomRows, n);
    yq = [topYPad; yq; bottomYPad];
end

function [params, values] = parseVarargin(args, defaults)
    % Parse variable arguments into parameter-value pairs.
    params = defaults(:, 1); values = defaults(:, 2);
    % (Implementation: Loop through args to match and assign.)
end