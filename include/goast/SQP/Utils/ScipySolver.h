
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <goast/Core.h>
#include <iostream>

namespace py = pybind11;

/*

// Convert Eigen matrix to NumPy array
py::array_t<double> eigen_matrix_to_numpy(const Eigen::MatrixXd &mat) {
    return py::array_t<double>({mat.rows(), mat.cols()}, mat.data());
}

// Convert NumPy array back to Eigen matrix
Eigen::MatrixXd numpy_to_eigen_matrix(const py::array_t<double> &array) {
    py::buffer_info buf_info = array.request();
    double *ptr = static_cast<double *>(buf_info.ptr);
    Eigen::Map<Eigen::MatrixXd> mat(ptr, buf_info.shape[0], buf_info.shape[1]);
    return mat;
}

// Function to solve a linear system using SciPy
Eigen::VectorXd solve_with_scipy(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Convert Eigen matrix and vector to NumPy arrays
    py::array_t<double> A_np = eigen_matrix_to_numpy(A);
    py::array_t<double> b_np = eigen_matrix_to_numpy(b);

    // Import scipy.linalg.solve in Python
    py::module scipy_linalg = py::module::import("scipy.linalg");

    // Solve the system using scipy.linalg.solve
    py::object result = scipy_linalg.attr("solve")(A_np, b_np);

    // Convert the result back to Eigen::VectorXd
    Eigen::VectorXd solution = numpy_to_eigen_matrix(result).col(0);
    return solution;
}
*/

class ScipySolver{
    private:
        py::module np;
        py::module scipy_linalg;

        std::vector<double> buffer_A;
        std::vector<double> buffer_rhs;

        // Convert Eigen matrix to NumPy array (pre-allocated memory)
        py::array_t<double> eigen_matrix_to_numpy(const Eigen::MatrixXd &mat, std::vector<double>& buffer) {
            std::memcpy(buffer.data(), mat.data(), mat.size() * sizeof(double));  // Copy data to buffer
            return py::array_t<double>({mat.rows(), mat.cols()}, buffer.data());  // Pass buffer data
        }

        // Convert NumPy array back to Eigen matrix
        Eigen::MatrixXd numpy_to_eigen_matrix(const py::array_t<double> &array) {
            py::buffer_info buf_info = array.request();
            double *ptr = static_cast<double *>(buf_info.ptr);
            Eigen::Map<Eigen::MatrixXd> mat(ptr, buf_info.shape[0], buf_info.shape[1]);
            return mat;
        }

    public:
        ScipySolver(size_t n){
            Py_Initialize();
            np = py::module::import("numpy");
            scipy_linalg = py::module::import("scipy.linalg");

            buffer_A.reserve(n*n);
            buffer_rhs.reserve(n);
        }

        ~ScipySolver()
        {
            buffer_A.clear();
            buffer_A.shrink_to_fit();
            buffer_rhs.clear();
            buffer_rhs.shrink_to_fit();
            Py_Finalize();
        }
        
        // Function to solve a linear system using SciPy
        Eigen::VectorXd solve_with_scipy(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
            // Convert Eigen matrix and vector to NumPy arrays
            py::array_t<double> A_np = eigen_matrix_to_numpy(A, buffer_A);
            py::array_t<double> b_np = eigen_matrix_to_numpy(b, buffer_rhs);

            /*
            py::tuple shape = A_np.attr("shape");
                
            // Get the number of rows (and columns for a square matrix)
            int n = shape[0].cast<int>();  // Number of rows
            int m = shape[1].cast<int>();  // Number of columns (use this for non-square matrices)
                
            // Create the identity matrix of appropriate size
            py::object identity_matrix = np.attr("eye")(n, m);  // Use n, m for general case
            py::object regularizer_obj = py::float_(1e-6);

            while(true)
            {
                py::object cond_value = np.attr("linalg").attr("cond")(A_np);
                //rcond = rcond_value.cast<double>();
                if(cond_value <= py::float_(1e-12))
                {
                    break;
                }
                // Add regularizer * Identity to A_128
                A_np = A_np + regularizer_obj * identity_matrix;
                regularizer_obj*=py::float_(10.0);
            }
                */

            // Solve the system using scipy.linalg.solve
            py::object result = scipy_linalg.attr("solve")(A_np, b_np);

            // Convert result back to double precision before returning
            py::object result_double = np.attr("array")(result, "float64");

            // Convert back to Eigen::VectorXd
            Eigen::VectorXd solution = numpy_to_eigen_matrix(result_double).col(0);
            return solution;
        }
};
