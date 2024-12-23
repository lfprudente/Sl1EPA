# **A Smooth ℓ₁-Exact Penalty Algorithm for Riemannian Optimization Problems**

This repository contains MATLAB codes related to the numerical experiments from the article:

> **E. G. Birgin, O. P. Ferreira, G. Haeser, N. Maculan, L. M. Ramirez, and L. F. Prudente,**  
> *Smoothing ℓ₁-Exact Penalty Method for Intrinsically Constrained Riemannian Optimization Problems*, 2024.

**Date:** December 2024

---

## 📂 **Contents**

**Main Routine**  
- **`main.m`**  
  - Controls the execution of numerical experiments.  
  - Allows users to select:  
    - The optimization problem to be solved.  
    - The smoothing function to be considered.

**Algorithm**  
- **`Sl1EPA.m`**  
  - Implements the **Smooth ℓ₁-Exact Penalty Algorithm**.

**Problem Routines**  
- **`evalf.m`**: Implements the **objective function**.  
- **`evalg.m`**: Implements the **Euclidean gradient** of the objective function.  
- **`evalcc.m`**: Implements the **constraints**.  
- **`evalnc.m`**: Implements the **Euclidean gradients** of the constraints.  

**Smoothing Routines**  
- **`phi.m`**: Implements the **smoothing functions**.  
- **`dphi.m`**: Implements the **derivatives** of the smoothing functions.  

---

## 🛠️ **Instructions**

### **Installing Dependencies**  
Make sure you have all required dependencies installed. Run the following command in the main folder:  

```sh
install_dependencies
```

### **Running the Numerical Experiments**  
To run the experiments:  
1. Open the file **`main.m`**.  
2. Select the **desired problem** and **smoothing function** in the script.  
3. Run the script and observe the results displayed in the **MATLAB command window**.

---

## 📦 **Third-Party Codes**

This repository includes third-party free software:

1. **Manopt 7.1**  
    - **Description:** MATLAB toolbox for optimization on manifolds.  
    - **Author:** Nicolas Boumal  
    - **Website:** [https://www.manopt.org/](https://www.manopt.org/)  
    - **License:** GNU General Public License (GPL) version 3  

⚠️ **Note:** Please review the specific licenses of the third-party codes before modifying or redistributing them. License information can be found in code comments or separate text files in their corresponding directories.

---

## 📄 **License**

This project is licensed under the **GNU General Public License (GPL) version 3**. For more details, see the **LICENSE** file included in this repository.

---

## 📚 **References**

- Boumal, N. *Manopt: A MATLAB Toolbox for Optimization on Manifolds*. [Manopt Documentation](https://www.manopt.org/).
