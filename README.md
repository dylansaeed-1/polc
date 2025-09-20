# POLC: Physics-in-the-Loop Learned Correction for 2D Darcy Flow

This project explores **Physics-in-the-Loop Learned Correction (POLC)** — a hybrid approach to surrogate modeling for PDEs. We focus on steady, incompressible Darcy flow in 2D porous media.

---

## What we’re doing
We train fast neural surrogates (U-Net, FNO) to map permeability fields **K(x,y)** → pressure **p(x,y)**.  
Instead of relying purely on data or physics losses, we insert **one lightweight correction step** during training and inference:

1. Predict pressure with a neural operator.  
2. Query a C++ finite-volume/finite-element solver (via gRPC) for the discrete residual and one conjugate-gradient update.  
3. Apply a learnable “trust” correction:  
   $
   p_\theta^{corr} = p_\theta + \alpha \Delta p
   $

---

## Why it matters
- **Conservation**: Better mass balance and boundary flux consistency.  
- **Robustness**: Improved generalization to out-of-distribution permeability fields.  
- **Efficiency**: Near the cost of a surrogate, far cheaper than iterative solvers.  
- **Modern angle**: Embeds numerical physics into ML, bridging PINNs and classical solvers.

---

## Approach
- **Data**: Synthetic permeability fields (Gaussian random fields, channelized facies, binary inclusions).  
- **Models**: U-Net and Fourier Neural Operator baselines.  
- **Solver coupling**: Residuals and optional CG updates via C++ backend.  
- **Uncertainty quantification**: Deep ensembles and residual-aware conformal calibration.  
- **Evaluation**: Accuracy, conservation, OOD robustness, runtime overhead, and calibration quality.

---

## Quickstart
```bash
# Setup environment
conda env create -f env.yml && conda activate polc

# Build C++ solver
cd cpp_solver && cmake -S . -B build && cmake --build build -j

# Train baseline
python ml/train_supervised.py --model unet --data data/grf_128

# Train with POLC (solver must be running)
python ml/train_polc.py --model unet --grpc localhost:50051
