% Load or generate correlation matrix
correlationMatrix = observedSignificanceMatrix;
% Set threshold for correlation-based network
%correlationThreshold = 0.5; % Adjust as needed
% Threshold the correlation matrix
%correlationMatrixThresholded = correlationMatrix > correlationThreshold;
% Create graph object for correlation-based network
correlationGraph = graph(correlationMatrix);
% Compute degree distribution
degreeDistribution = degree(correlationGraph);
% Calculate degree distribution probability P(k)
maxDegree = max(degreeDistribution);
degreeCounts = histcounts(degreeDistribution, 0:maxDegree);
degreeProbabilities = degreeCounts / sum(degreeCounts);
% Plot degree distribution probability
figure;
plotted_bar = bar(0:maxDegree - 1, degreeProbabilities);
xlabel('Degree (k)');
ylabel('Probability P(k)');
title('Degree Distribution Probability P(k) of Correlation-based Network');
X = plotted_bar.XData;
X = transpose(X);
Y = plotted_bar.YData;
Y = transpose(Y);

% Perform a power-law fit

%fo = fitoptions('Method','NonlinearLeastSquares',...
 %              'StartPoint',[1 -1]);

fo = fitoptions('Method','NonlinearLeastSquares');

%model = fittype('a * x^b', 'independent', 'x', 'dependent', 'y', 'options', fo);
model = fittype('C * x^(-gamma)', 'independent', 'x', 'dependent', 'y', 'options', fo);

powerLawFit = fit(X, Y, model);

% Get the fit parameters (C and gamma)
C = powerLawFit.C;
gamma = powerLawFit.gamma;

% Assess the goodness of fit
% For example, calculate R-squared
rsquared = powerLawFit.rsquare;
fprintf('R-squared: %f\n', rsquared);

% Plot the data and the power-law fit
scatter(X, Y, 'b', 'filled');
hold on;
plot(powerLawFit, 'r');
xlabel('X');
ylabel('Y');
title('Power Law Fit');
legend('Data', 'Power Law Fit');



