# Income and Wealth Distribution using Julia
This repo replicate the `diffussion` and `fat tail` income and wealth distribution in a Aiyagari-Bewley-Huggett model from [Income and Wealth Distribution in Macroeconomics: A Continuous-Time Approach](https://benjaminmoll.com/wp-content/uploads/2019/07/HACT.pdf). 

The Julia codes replicates in a efficient and intuitive way the steady state. These are an adaptation of the original codes in [Matlab](https://benjaminmoll.com/codes/). The repo is divided into two folders: **diffusion** and **fat_tail**. Each has the following structure: 

- `1_Parameters.jl`: definition of the parameters.
- `2_Steady_state.jl`: i) steady state, ii) partial equilibrium, iii) Fokker-Planck equation, v) updating of prices, and vi) individual solution of the value function.  
-  `3_Main.ipynb`: model solution.

