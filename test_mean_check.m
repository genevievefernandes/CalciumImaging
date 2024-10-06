% Sample matrix, replace this with your actual data
matrix = randi([1, 10], 10, 2);

% Set the interval size
interval_size = 2;

% Get the number of rows and columns in the matrix
num_rows = size(matrix, 1);
num_columns = size(matrix, 2);

% Initialize a vector to store the means
means = zeros(ceil(num_rows / interval_size), num_columns);

% Iterate over rows and calculate mean for each interval
for i = 1:interval_size:num_rows
    end_index = min(i + interval_size - 1, num_rows);
    interval_data = matrix(i:end_index, :);
    interval_mean = mean(interval_data, 1);
    means(ceil(i / interval_size), :) = interval_mean;
end

% Now 'means' contains the mean values for each interval
% Each row in 'means' corresponds to the mean of the corresponding interval

prepared_matrix = 1 ./ means;