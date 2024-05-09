# Multi-Cost-SVM (and Probabilistic Safety Regions for exponential distributions)

Multi Cost SVM (MC-SVM) is a variant of Support Vector Machines (SVM) designed to accommodate multiple cost scenarios. By introducing multiple weighting parameters $\tau$,  MC-SVM adapts the cost function to balance false positive and false negative errors, enhancing the model's robustness across diverse scenarios. The result is a separation hyperplane indipendent from the sample probability of the data.

$$
\begin{align}
& &\underset{\mathbf{w}, \mathbf{c}, \boldsymbol{\xi}_1, \ldots, \boldsymbol{\xi}_n}{\text{minimize}} \frac{1}{2\eta} 
\mathbf{w}^\top \mathbf{w} + \frac{1}{2} 
\end{align}
$$

### Key Features:
_Parameterized Cost Function_: MC-SVM incorporates a parameter $\tau$ to influence the cost function's behavior towards different types of errors. This parameterization allows to weight the SVMs with different weighting parameters, reducing the unbalanceness of the data and helping training a more robust algorithm.
_System of SVMs_: The algorithm constructs a system of $m$ SVMs using the same dataset but varying weights and offsets. Each SVM corresponds to a different value of $\tau$, enabling the model to adapt to various cost scenarios.

The optimization problem is solved in its dual form

Usage:
To utilize MC-SVM in your projects, follow these steps:

Download the Code: Clone the repository containing the MC-SVM implementation.
Configure Parameters: Adjust the value of 
𝜏
τ and other parameters according to your application requirements.
Train the Model: Provide your dataset and train the MC-SVM model using the provided training routines.
Evaluate Performance: Evaluate the model's performance on your test dataset and analyze its behavior under different cost scenarios.
Example:
python
Copy code
# Import MC-SVM module
from mc_svm import MCSVM

# Initialize MC-SVM with desired parameters
mc_svm = MCSVM(tau=0.5, eta=1.0)

# Train the model on the training dataset
mc_svm.train(X_train, y_train)

# Evaluate the model on the test dataset
accuracy = mc_svm.evaluate(X_test, y_test)
print("Accuracy:", accuracy)
Visualization:
Visualize the algorithm's performance across different 
𝜏
τ values using provided plotting utilities. These visualizations can offer insights into the algorithm's behavior and guide parameter selection.


Contributions and Feedback:
Contributions to the MC-SVM algorithm are welcome! Feel free to submit bug reports, feature requests, or pull requests to improve the algorithm's functionality and usability.

References:
Koenker, R. (2005). Quantile regression. Cambridge university press.
Note: This README provides a brief overview of the MC-SVM algorithm and its usage. For detailed instructions and documentation, refer to the provided codebase and accompanying documentation files.
