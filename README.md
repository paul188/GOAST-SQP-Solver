# Optimization of Folds and Developability on Discrete Plates

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Thesis](https://img.shields.io/badge/Thesis-Mathematics-green)](https://www.ins.uni-bonn.de/)

This repository contains the implementation of a specialized **Sequential Quadratic Programming (SQP) Solver** for the optimization of folds and developability on discrete plates. This work was developed as part of a Master's Thesis in Mathematics at the **University of Bonn**.

## 📖 Project Overview

This project extends the [GOAST (Geodesic Orbitals and Small Transformations)](https://gitlab.com/numod/goast) library to handle complex geometric optimization problems. The core focus is on:

* **Developable Surfaces:** Ensuring that discrete meshes can be flattened into a plane without stretching or tearing.
* **Fold Optimization:** Strategically placing and adjusting folds to reach target 3D geometries.
* **SQP Solver:** A custom Sequential Quadratic Programming implementation designed to handle high-dimensional non-linear constraints inherent in discrete differential geometry.

---

## 🛠 Features

* **Custom SQP Implementation:** Tailored for energy functionals related to discrete plate theory.
* **Constraint Management:** Efficient handling of isometric and developability constraints.
* **GOAST Integration:** Built directly upon the robust GOAST framework for geodesic calculations and manifold optimization.
* **Visualization Support:** Tools to export optimization steps for analysis in software like ParaView or Polyscope.

---

## 🖼 Gallery

> **Note:** To display images, save your thesis plots as `.png` files in a folder named `docs/` and update the links below.

| Initial Mesh | SQP Optimization Progress | Final Folded State |
| :---: | :---: | :---: |
| ![Initial State](docs/initial.png) | ![Process](docs/process.gif) | ![Final State](docs/final.png) |

---

## ⚙️ Installation & Build

### Prerequisites
* C++17 compatible compiler (GCC, Clang, or MSVC)
* [CMake](https://cmake.org/) (>= 3.14)
* [Eigen3](https://eigen.tuxfamily.org/)
  
## 🚀 Quick Start

To run a fold optimization experiment using the SQP solver:

```bash
./bin/fold_optimizer --config configs/your_experiment.json ```bash

## 📊 Results

The SQP solver implemented in this repository has been tested on a variety of discrete plate optimization problems. Below are the key results highlighting the solver's performance and capabilities.

### Fold Optimization on Discrete Plates
The solver successfully identifies optimal fold configurations that minimize energy functionals while strictly adhering to developability and boundary constraints.

| Initial State | SQP Convergence | Final Geometry |
| :---: | :---: | :---: |
| ![Flat Plate](docs/img/initial_plate.png) | ![Energy Decay](docs/img/convergence_plot.png) | ![Folded Result](docs/img/optimized_fold.png) |

### Key Performance Metrics
* **Convergence Efficiency**: The SQP approach demonstrates superior convergence rates compared to standard gradient-based methods for non-linear developability constraints.
* **Constraint Satisfaction**: Achieves high-precision isometric mapping, ensuring the discrete plates remain developable throughout the optimization process.
* **Versatility**: Capable of handling complex target curvatures and intricate fold patterns on high-resolution meshes.

### Example Cases
* **Target Curvature Matching**: Optimizing a flat sheet to match a specified Gaussian curvature distribution while maintaining a discrete plate structure.
* **Boundary-Driven Folds**: Generating folding patterns based on external boundary forces and prescribed fold line locations.
