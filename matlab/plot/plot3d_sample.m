% Define the range of x and y
% meshgrid generates a grid of x and y values
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);

% Define the function z = sin(sqrt(x^2 + y^2))
% '.' before '^' and '*' makes the operations element-wise
z = sin(sqrt(x.^2 + y.^2));

% Create the 3D plot
% surf generates the 3D surface plot
surf(x, y, z)