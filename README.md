This file can perform some simple mathematical calculation operations, such as basic operations, data statistics, matrix operations, prime number operations, as well as array operations
File 1: mathlib_matrix. c (matrix calculation library)
effect:
Provide basic operations and calculation functions for matrices, including matrix creation, destruction, addition, subtraction, multiplication, transposition, determinant calculation, matrix power, inverse matrix, adjoint matrix, matrix rank, eigenvalue calculation, etc.
Main functions:
create_matrix： Create a matrix with a specified number of rows and columns.
destroy_matrix： Destroy the matrix and free up the memory it occupies.
matrix_add： Matrix addition.
matrix_subtract： Matrix subtraction.
matrix_multiply： Matrix multiplication.
matrix_transpose： Matrix transpose.
matrix_determinant： Calculate the determinant of a matrix (only supports 2x2 and 3x3 matrices).
matrix_trace： Calculate the trace of the matrix.
matrix_power： Calculate the power of a matrix.
matrix_inverse： Compute the inverse of a matrix (only supports 2x2 and 3x3 matrices).
matrix_adjoint： Calculate the adjoint matrix of the matrix (only supports 2x2 and 3x3 matrices).
matrix_rank： Calculate the rank of the matrix (only supports 2x2 and 3x3 matrices).
matrix_eigenvalues： Calculate the eigenvalues of the matrix (only supports 2x2 and 3x3 matrices).
File 2: mathlib_prime. c (library of prime numbers and sorting algorithms)
effect:
Provide calculation functions related to prime numbers (such as detecting prime numbers, counting prime numbers, filtering prime number intervals) and the implementation of various sorting algorithms.
Main functions:
is_prime： Check if a number is prime.
prime_count： Calculate the number of prime numbers that are less than or equal to a certain number.
prime_sieve： Filter out all prime numbers within the specified interval.
quick_sort： Quick sort algorithm.
merge_sort： Merge sort algorithm.
heap_sort： Heap sorting algorithm.
radix_sort： Cardinality sorting algorithm.
File 3: mathlib_stat. c (Statistical Computing Library)
effect:
Provide statistical calculation functions, including mean, variance, standard deviation, covariance, correlation coefficient, skewness, kurtosis, median, mode, range, quartiles, normal distribution CDF, z-score, binomial probability, Poisson probability, geometric probability, hypergeometric probability, linear regression coefficient, coefficient of determination (R-squared), mean absolute error, mean square error, root mean square error, mean absolute percentage error, etc.
Main functions:
calculate_mean： Calculate the mean.
calculate_variance： Calculate variance.
calculate_standard_deviation： Calculate the standard deviation.
calculate_covariance： Calculate covariance.
calculate_correlation： Calculate the correlation coefficient.
calculate_skewness： Calculate skewness.
calculate_kurtosis： Calculate kurtosis.
calculate_median： Calculate the median.
calculate_mode： Calculate the mode.
calculate_range： Calculate the range.
calculate_quartiles： Calculate the quartiles.
calculate_normal_cdf： Calculate the cumulative distribution function (CDF) of a normal distribution.
calculate_z_score： Calculate the z-score.
calculate_binomial_probability： Calculate binomial probability.
calculate_poisson_probability： Calculate Poisson probability.
calculate_geometric_probability： Calculate geometric probability.
calculate_hypergeometric_probability： Calculate the hypergeometric probability.
calculate_linear_regression： Calculate the linear regression coefficients.
calculate_r_squared： Calculate the coefficient of determination (R-squared).
calculate_mae： Calculate the average absolute error.
calculate_mse： Calculate the mean square error.
calculate_rmse： Calculate the root mean square error.
calculate_mape： Calculate the average absolute percentage error.
File 4: mathlib array. c (dynamic array operation library)
effect:
Implement a dynamic array structure and provide various operations on the dynamic array, such as adding, inserting, deleting, retrieving, setting elements, searching, sorting, etc.
Main functions:
create_dynamic_array： Create a dynamic array.
destroy_dynamic_array： Destroy dynamic arrays.
dynamic_array_add： Add elements to a dynamic array.
dynamic_array_insert： Insert elements at the specified position in a dynamic array.
dynamic_array_remove： Delete the element at the specified position in the dynamic array.
dynamic_array_get： Retrieve the element at the specified position in a dynamic array.
dynamic_array_set： Set the element at the specified position in a dynamic array.
dynamic_array_size： Get the size of the dynamic array.
dynamic_array_capacity： Retrieve the capacity of a dynamic array.
dynamic_array_find： Find the position of an element in a dynamic array.
dynamic_array_bubble_sort： Bubble sort dynamic arrays.
dynamic_array_selection_sort： Sort and select dynamic arrays.
dynamic_array_insertion_sort： Insert sort dynamic arrays.
dynamic_array_shell_sort： Perform Hill sort on dynamic arrays.
dynamic_array_merge_sort： Merge and sort dynamic arrays.
dynamic_array_heap_sort： Sort dynamic arrays by heap.
dynamic_array_radix_sort： Sort dynamic arrays by cardinality.
File 5: mathlib. base. c (Basic Mathematical Computation Library)
effect:
Provide basic mathematical calculation functions, including factorial, greatest common divisor, least common multiple, Fibonacci sequence, combination number, permutation number, power operation, square root, logarithm, exponent, trigonometric function, hyperbolic function, inverse trigonometric function, inverse hyperbolic function, absolute value, rounding, square, cubic, and array versions of these functions.
Main functions:
factorial： Calculate factorial.
gcd： Calculate the greatest common divisor.
lcm： Calculate the least common multiple.
fibonacci： Calculate the Fibonacci sequence.
combination： Calculate the number of combinations.
permutation： Calculate the number of permutations.
power： Calculate power.
sqrt_custom： Calculate the square root.
ln： Calculate natural logarithm.
log10_custom： Calculate common logarithms.
exp_custom： Calculate the exponential function.
sin_custom： Calculate the sine function.
cos_custom： Calculate the cosine function.
tan_custom： Calculate the tangent function.
sinh_custom： Calculate hyperbolic sine function.
cosh_custom： Calculate the hyperbolic cosine function.
tanh_custom： Calculate the hyperbolic tangent function.
asin_custom： Calculate the inverse sine function.
acos_custom： Calculate the inverse cosine function.
atan_custom： Calculate the arctangent function.
asinh_custom： Calculate the inverse hyperbolic sine function.
acosh_custom： Calculate the inverse hyperbolic cosine function.
atanh_custom： Calculate the inverse hyperbolic tangent function.
abs_custom： Calculate the absolute value.
floor_custom： Calculate the rounding function.
ceil_custom： Calculate the top function.
round_custom： Calculate the rounding function.
square： Calculate the square function.
cube： Calculate the cubic function.
Array version functions, such as factorial'array, fibonacci_array, combination'array, etc., are used to generate an array of calculation results.
