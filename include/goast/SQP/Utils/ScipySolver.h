
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <goast/Core.h>
#include <iostream>

namespace py = pybind11;

class __attribute__((visibility("default"))) ScipySolver{
    public:
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

            //pybind11::scoped_interpreter guard{}; // Start the interpreter and keep it alive

            //setenv("PYTHONHOME", "/home/paul_johannssen/anaconda3/envs/myenv", 1);
            //setenv("PYTHONPATH", "/home/paul_johannssen/anaconda3/envs/myenv/lib/python3.13/site-packages", 1);
            
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
            np = py::module();
            scipy_linalg = py::module();
            Py_Finalize();
        }
        
        // Function to solve a linear system using SciPy
        Eigen::VectorXd solve_with_scipy(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
            // Convert Eigen matrix and vector to NumPy arrays
            py::array_t<double> A_np = eigen_matrix_to_numpy(A, buffer_A);
            py::array_t<double> b_np = eigen_matrix_to_numpy(b, buffer_rhs);

            // Solve the system using scipy.linalg.solve
            py::object result = scipy_linalg.attr("solve")(A_np, b_np);

            // Convert result back to double precision before returning
            py::object result_double = np.attr("array")(result, "float64");

            // Convert back to Eigen::VectorXd
            Eigen::VectorXd solution = numpy_to_eigen_matrix(result_double).col(0);
            return solution;
        }

        Eigen::VectorXd solve_with_scipy_pseudoinv(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
            // Convert Eigen matrix and vector to NumPy arrays
            py::array_t<double> A_np = eigen_matrix_to_numpy(A, buffer_A);
            py::array_t<double> b_np = eigen_matrix_to_numpy(b, buffer_rhs);

            // Solve the system using scipy.linalg.solve
            py::dict kwargs;
            kwargs["atol"] = 1e-6;
            py::object pseudo_inv = scipy_linalg.attr("pinv")(A_np, **kwargs);

            // Convert result back to double precision before returning
            py::object result_double = np.attr("array")(pseudo_inv, "float64");

            // Convert back to Eigen::VectorXd
            Eigen::MatrixXd pseudoinv = numpy_to_eigen_matrix(result_double);
            Eigen::VectorXd solution = pseudoinv.transpose() * b;
            return solution;
        }
};
