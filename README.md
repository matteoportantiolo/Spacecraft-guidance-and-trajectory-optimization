# Spacecraft Guidance and Trajectory Optimization

## Overview
This project addresses three trajectory and guidance problems in astrodynamics:

1. **Preliminary trajectory analysis**
2. **Orbit transfer via shooting methods**
3. **Low-thrust optimal control using Pontryagin Maximum Principle**

The work combines numerical optimization, orbital dynamics, and optimal control techniques to solve realistic spacecraft guidance problems.

---

## 1 – Preliminary Trajectory Analysis

A preliminary analysis is performed to determine feasible initial conditions and transfer timing:

- Identification of transfer windows based on orbital geometry
- Propagation of trajectories using simplified and n-body models
- Comparison of dynamical models to assess perturbation effects

This step provides consistent initial guesses and physical insight for the subsequent optimization problems.

---

## 2 – Trajectory Optimization via Shooting Methods

### Problem Description
A spacecraft transfer trajectory is designed between two states using high-fidelity dynamical models. The goal is to compute a feasible and optimal trajectory satisfying boundary conditions on position and velocity.

### Methods

#### Shooting Methods
The boundary value problem is solved using:

- **Single Shooting**
  - Direct root-finding on initial conditions
  - Fast but sensitive to initialization

- **Multiple Shooting**
  - Time domain split into segments
  - Improved robustness and convergence
  - Higher computational cost

#### Dynamical Models
- Restricted dynamical model (PBRFBP)
- High-fidelity **n-body propagation**

A comparison between models highlights the impact of perturbations on trajectory evolution.

### Results

- Multiple shooting significantly improves solution quality:
  - Lower ΔV
  - Shorter transfer time
- The optimal trajectory approaches a **Pareto-optimal solution**
- n-body propagation reveals non-planar effects and long-term divergence

---

## 3 – Low-Thrust Optimal Control (PMP)

### Problem Description
A low-thrust spacecraft performs an orbit raising maneuver:
- From 800 km circular orbit
- To 1000 km circular orbit with inclination change

The objective is to **minimize exposure to orbital debris**, modeled as a spatial density function.

### Optimal Control Formulation

The problem is solved using **Pontryagin Maximum Principle (PMP)**:

- State dynamics and costate equations derived
- Hamiltonian constructed
- Boundary conditions enforced via shooting

The optimal control structure is:

- **Thrust direction → aligned with primer vector**
- **Thrust magnitude → determined by switching function**

### Numerical Solution

- Nonlinear root-finding on:
  - Initial costates
  - Final time
- Random initialization for robustness
- Numerical continuation for varying thrust levels

### Results

- Optimal trajectories balance:
  - Transfer time
  - Fuel consumption
  - Debris exposure
- Lower thrust leads to:
  - Longer transfers
  - Different control structure

---

## Key Techniques and Contributions

- Boundary value problem solving via shooting methods
- Multiple shooting for improved convergence
- Optimal control using PMP
- Primer vector interpretation of optimal thrust
- Numerical continuation strategies
- Comparison between simplified and high-fidelity dynamics

---

## Implementation

The project includes:
- Nonlinear solvers for shooting methods
- Orbital propagation (2-body and n-body)
- Optimal control formulation and solution
- Trajectory visualization and analysis tools

---

## Key Concepts

- Trajectory optimization
- Shooting methods (single vs multiple)
- Pontryagin Maximum Principle
- Low-thrust guidance
- Primer vector theory
- Orbital dynamics and perturbations

---

## Author
Matteo Portantiolo  
MSc Space Engineering – GNC
