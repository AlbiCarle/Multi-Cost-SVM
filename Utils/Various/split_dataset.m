function [X_train, Y_train, X_test, Y_test, X_calib, Y_calib] = split_dataset(X, Y, train_size, test_size, calib_size)
    % Input:
    % X: Input features matrix (each row is a sample, each column is a feature)
    % Y: Target values vector (column vector corresponding to the samples)
    % train_size: Size of the training set (number of samples for training)
    % test_size: Size of the test set (number of samples for testing)
    % calib_size: Size of the calibration set (number of samples for calibration)
    %
    % Output:
    % X_train, Y_train: Training set (features and target values)
    % X_test, Y_test: Test set (features and target values)
    % X_calib, Y_calib: Calibration set (features and target values)

    % Check if the specified sizes are valid
    if train_size <= 0 || test_size <= 0 || calib_size <= 0
        error('Sizes of training, test, and calibration sets must be positive integers.');
    end
    
    total_samples = size(X, 1);
    total_set_size = train_size + test_size + calib_size;

    % Check if the specified sizes exceed the total number of samples
    if total_set_size > total_samples
        error('Sizes of training, test, and calibration sets exceed the total number of samples.');
    end

    % Randomly permute the indices of the samples
    perm_indices = randperm(total_samples);

    % Split the dataset based on the specified sizes
    X_train = X(perm_indices(1:train_size), :);
    Y_train = Y(perm_indices(1:train_size), :);

    X_test = X(perm_indices(train_size+1:train_size+test_size), :);
    Y_test = Y(perm_indices(train_size+1:train_size+test_size), :);

    X_calib = X(perm_indices(train_size+test_size+1:train_size+test_size+calib_size), :);
    Y_calib = Y(perm_indices(train_size+test_size+1:train_size+test_size+calib_size), :);
end
