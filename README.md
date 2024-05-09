# Multi-Cost-SVM (and Probabilistic Safety Regions for exponential distributions)

Multi Cost SVM (MC-SVM) is a variant of Support Vector Machines (SVM) designed to accommodate multiple cost scenarios. By introducing multiple weighting parameters $\tau$,  MC-SVM adapts the cost function to balance false positive and false negative errors, enhancing the model's robustness across diverse scenarios. The result is a separation hyperplane indipendent from the sample probability of the data.

### Key Features:
_Parameterized Cost Function_: MC-SVM incorporates a parameter $\tau$ to influence the cost function's behavior towards different types of errors. This parameterization allows to weight the SVMs with different weighting parameters, reducing the unbalanceness of the data and helping training a more robust algorithm.
_System of SVMs_: The algorithm constructs a system of $m$ SVMs using the same dataset but varying weights and offsets. Each SVM corresponds to a different value of $\tau$, enabling the model to adapt to various cost scenarios.

The optimization problem is solved in its dual form

Usage:
To utilize MC-SVM in your projects, follow these steps:

_Download the Code_: Clone the repository containing the MC-SVM implementation.
_Configure Parameters_: Adjust the value of $\tau$, the kernels and other parameters according to your application requirements.
_Train the Model_: Provide your dataset and train the MC-SVM model using the provided training algorithm.
_Evaluate Performance_: Evaluate the model's performance on your test dataset and analyze its behavior under different cost scenarios.

### Example:
Matlab

_For a dataset composed by data sampled with different probabilities_

Tau  = rand(1,9); 
m = size(Tau,2);

kernel = 'polynomial';
param = 3;
eta = .001;

alpha_bar = MCSVM_Train(Xtr, Ytr, kernel, param, Tau, eta); # best hyperplane common to all the data

_Specializing to a dataset with a known (or estimated) sample probability_

tau = 1-epsilon; # to control the false positives 

alpha_c = SSVM_Train_c(Xtr, Ytr, Xcl_p, Ycl_p, kernel, param, tau, eta, alpha_bar);

b = offset_c(Xtr, Ytr, Xcl_p, Ycl_p, alpha_c, kernel, param, eta, tau, alpha_bar); # best offset that realizes the control of the false positive ration on the desired (calibration) set.

_test_

y_pred_ts = SSVM_Test(Xtr, Ytr, Xts_p, alpha_bar, b, 0, kernel, param, eta);

[TPR_SSVM, FPR_SSVM, TNR_SSVM, FNR_SSVM, F1_SSVM, ACC_SSVM] = ConfusionMatrix(Yts_p, y_pred_ts,'on');

disp(['False positive rate:',num2str(FPR_SSVM)])

Contributions and Feedback:
Contributions to the MC-SVM algorithm are welcome! Feel free to submit bug reports, feature requests, or pull requests to improve the algorithm's functionality and usability.

References:

