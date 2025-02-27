#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Dense>
#include <iostream>

namespace py = pybind11;

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

int main()
{
    Py_Initialize();
    if (!Py_IsInitialized()) {
        std::cerr << "Python initialization failed!" << std::endl;
        return 1;
    }
    PyObject *obj;
    Eigen::MatrixXd A(2, 2);
    A << 3, 1,
        1, 2;

    Eigen::VectorXd b(2);
    b << 9, 8;

    // Call the C++ function (which calls SciPy in Python)
    Eigen::VectorXd solution = solve_with_scipy(A, b);
    std::cout<<"Solution: "<<solution[0]<<"; "<<solution[1]<<std::endl;
    return 0;
}