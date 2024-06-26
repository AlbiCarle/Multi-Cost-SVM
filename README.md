# Multi-Cost-SVM (and Probabilistic Safety Regions for exponential distributions)

<img src = Images/coolGaussians.png width ="800">

Multi Cost SVM (MC-SVM) is a variant of Support Vector Machines (SVM) designed to accommodate multiple cost scenarios. By introducing multiple weighting parameters $\tau$,  MC-SVM adapts the cost function to balance false positive and false negative errors, enhancing the model's robustness across diverse scenarios. The result is a separation hyperplane indipendent from the sample probability of the data.

<img src = Images/minimum1.png width="300">

This algorithm was inspired by the concept of _Probabilistic Safety Region_ (PSR) 

<div style="text-align:center;">
    <img src="Images/Phi.png" width="200">
</div>

i.e., the region where in high probability is possible to observe the event $S$, that, we can suppose, represents a "safe" situation. It is interesting to note, and these considerations are reported in the code, that for exponential distributions the PSR takes the interesting form of a radius controllable set:

<div style="text-align:center;">
    <img src="Images/Phi2.png" width="200">
</div>

### Key Features:
__Parameterized Cost Function__: MC-SVM incorporates a parameter $\tau$ to influence the cost function's behavior towards different types of errors. This parameterization allows to weight the SVMs with different weighting parameters, reducing the unbalanceness of the data and helping training a more robust algorithm.

__System of SVMs__: The algorithm constructs a system of $m$ SVMs using the same dataset but varying weights and offsets. Each SVM corresponds to a different value of $\tau$, enabling the model to adapt to various cost scenarios.

The optimization problem is solved in its dual form

<img src = Images/minimum2.png width="300">

leading to the separation hyperplane 

<img src = Images/sep.png width="150">

The error in the prediction (false or negative ratio) is then controlled using the following algorithm, based on the quantile regression idea that, discarding the regularization parameter (possible because we computed an independent hyperplane with the algorithm above), the weighting parameter corresponds to the false negative ratio:

<img src = Images/Minimum3.png width="300">

### Usage:
To utilize MC-SVM in your projects, follow these steps:

__Download the Code__: Clone the repository containing the MC-SVM implementation.

__Configure Parameters__: Adjust the value of $\tau$, the kernels and other parameters according to your application requirements.

__Train the Model__: Provide your dataset and train the MC-SVM model using the provided training algorithm.

__Evaluate Performance__: Evaluate the model's performance on your test dataset and analyze its behavior under different cost scenarios.

### Example:

Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

_For a dataset composed by data sampled with different probabilities_

Tau  = rand(1,9); 
m = size(Tau,2);

kernel = 'polynomial';

param = 3;

eta = .001;

alpha_bar = MCSVM_Train(Xtr, Ytr, kernel, param, Tau, eta); # best hyperplane common to all the data

_Specializing to a dataset with a known (or estimated) sample probability_

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 1-epsilon; # to control the false positives 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_c = SSVM_Train_c(Xtr, Ytr, Xcl_p, Ycl_p, kernel, param, tau, eta, alpha_bar);

b = offset_c(Xtr, Ytr, Xcl_p, Ycl_p, alpha_c, kernel, param, eta, tau, alpha_bar); # best offset that realizes the control of the false positive ration on the desired (calibration) set.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

_test_

y_pred_ts = SSVM_Test(Xtr, Ytr, Xts_p, alpha_bar, b, 0, kernel, param, eta);

[TPR_SSVM, FPR_SSVM, TNR_SSVM, FNR_SSVM, F1_SSVM, ACC_SSVM] = ConfusionMatrix(Yts_p, y_pred_ts,'on');

disp(['False positive rate:',num2str(FPR_SSVM)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Contributions and Feedback:
Contributions to the MC-SVM algorithm are welcome! Feel free to submit bug reports, feature requests, or pull requests to improve the algorithm's functionality and usability.

### References:



