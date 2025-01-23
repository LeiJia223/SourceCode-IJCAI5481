# Operating Environment

Please ensure that the latest version of MATLAB is installed and includes all necessary toolboxes. All experiments were carried out in the MATLAB R2024a environment. The hardware platform is an Intel Core i7 - 11700 CPU with 32 GB of memory, and the operating system is Microsoft Windows 11.

# Directory Structure

- **EX1**: Contains the dataset and results of Experiment 1 in the paper (corresponding to Figure 1), including the source code for the comparison of different parameters in low-dimensional settings:
  - Group 1: γ = 1, δ = 0.5, α = 1
  - Group 2: γ = 10, δ = 0.5, α = 1
  - Group 3: γ = 10, δ = 0.9, α = 1
  - Group 4: γ = 10, δ = 0.9, α = 10

- **EX2**: Contains the dataset and results of Experiment 2 in the paper (corresponding to Figure 2), including the source code for the comparison of different parameters in high-dimensional settings:
  - Group 1: γ = 1, δ = 0.5, α = 1
  - Group 2: γ = 10, δ = 0.5, α = 1
  - Group 3: γ = 10, δ = 0.9, α = 1
  - Group 4: γ = 10, δ = 0.9, α = 10

- **EX3**: Contains the dataset and results of Experiment 3 in the paper (corresponding to Figure 3), including the source code for adding different types of noise to the model:
  - No noise
  - Constant noise
  - Time-varying noise
  - Random noise

- **EX4**: Contains the dataset and results of Experiment 4 in the paper (corresponding to Figure 4), including the source code for comparing five models under different low-dimensional noises:
  - No noise
  - Constant noise
  - Time-varying noise
  - Random noise

- **EX5**: Contains the dataset and results of Experiment 5 in the paper (corresponding to Figure 5), including the source code for comparing five models under different high-dimensional noises:
  - No noise
  - Constant noise
  - Time-varying noise
  - Random noise

- **Results**: Contains the `.fig` files of the output results.

# Experiment Details

## Experiment 1

All files of Experiment 1 are stored in the `EX1` folder. You can run the programs to verify the experimental results.

### Example Usage

Suppose you need to verify the contents in the `EX1` folder, you can run the program `com_para_2d.m` in the folder.

After the operation is completed, MATLAB will generate a figure, which shows the corresponding operation results.

## Experiment 2

All files of Experiment 2 are stored in the `EX2` folder. You can run the programs to verify the experimental results.

### Example Usage

Suppose you need to verify the contents in the `EX2` folder, you can run the program `com_para_10d.m` within the folder.

After the run is completed, MATLAB will generate a figure to display the corresponding operation results.

## Experiment 3

All files of Experiment 3 are stored in the `EX3` folder. You can run the programs to verify the experimental results.

### Example Usage

Suppose you need to verify the content in the `EX3` folder. For example, to verify the convergence of the NORNN model when there is no noise, you can uncomment the following line of code in the program `NORNN.m`:

```matlab
%Nt=0;

```
After the operation is completed, MATLAB will generate a figure to display the corresponding operation results.
## Experiment 4

All files of Experiment 4 are stored in the `EX4` folder. You can run the programs to verify the experimental results.

### Example Usage

To verify the convergence of five models under low-dimensional and noise-free conditions, uncomment the following line of code in the program `com_model_2d.m`:

```matlab
%Nt=0;
%dXdt_matrix_nornn = -1*n3 * Et - Pt_pt - dAt * X + Nt ; 
```
If you need to verify the convergence of the five models under low - dimensional time - varying noise, you can uncomment the following lines of code in the program `com_model_2d.m`:
```matlab
%Nt=0.1*sin(t);
%dXdt_matrix_nornn = -1*n3 * AFMSbp(Et) - Pt_pt - dAt * X + Nt 
```
After running, MATLAB will generate a figure to display the corresponding results.
## Experiment 5
All files of Experiment 5 are stored in the `EX5` folder. You can run the programs to verify the experimental results.
### Example Usage
Suppose you need to verify the content in the EX5 folder. For instance, to verify the convergence of the five models under high - dimensional and noise - free conditions, you can uncomment the following lines of code in the program `com_model_10d.m`:
```matlab
%Nt=0;
%dXdt_matrix_nornn = -1*n3 * Et - Pt_pt - dAt * X + Nt ;
```
If you need to verify the convergence of the five models under high - dimensional time - varying noise, you can uncomment the following lines of code in the program `com_model_10d.m`:
```matlab
%Nt=0.1*sin(t);
%dXdt_matrix_nornn = -1*n3 * AFMSbp(Et) - Pt_pt - dAt * X + Nt ;
```
After the run is completed, MATLAB will generate a figure to display the corresponding operation results.














