% Step 1: Create a CSV file with two columns of data
% This is done manually or programmatically. For example, the file 'data.csv' might contain:
% 1,2
% 3,4
% 5,6
% 7,8

% Step 2: Load the data from the CSV file
data = csvread('data.csv');

% Step 3: Store the first column of data in one variable and the second column in another variable
x = data(:, 1);
y = data(:, 2);

% Step 4: Plot the two-dimensional data
plot(x, y);
xlabel('x');
ylabel('y');
title('Plot of CSV Data');