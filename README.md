# Data Analysis Histograms Overview

This document provides a detailed explanation of the histograms used in the data analysis of **pp LHC23zzh pass4 (LHC23g3) MC anchored data**, including their classification into **Reconstructed Histograms** and **Simulated (MC Truth) Histograms**.

## Dataset Information
Analysis performed on:
- **AO2D file** from 2024 MC production
- Anchored to pp LHC23zzh pass4 (LHC23g3) dataset
- Uses both reconstructed tracks and MC truth information

## Reconstructed Histograms:
These histograms correspond to the **reconstructed tracks** from the detector, filled with track-level information:

- **eventCounter**: Number of processed events
- **etaHistogram**: η distribution of reconstructed tracks
- **ptHistogram**: Transverse momentum ($p_T$) distribution
- **ptResolution**: $p_T$ resolution (reconstructed vs generated)
- **tpcDedxVsPt**: TPC $dE/dx$ vs $p_T$ correlation
- **nSigmaPiVsPt**, **nSigmaKaVsPt**, **nSigmaPrVsPt**, **nSigmaElVsPt**: $n\sigma$ distributions for pions, kaons, protons, and electrons vs $p_T$
- **nSigmaPiVsP**: $n\sigma(\pi)$ vs $p$ distribution
- **ptSpectrumPion**, **ptSpectrumKaon**, **ptSpectrumProton**, **ptSpectrumElectron**: $p_T$ spectra using 3σ PID cuts
- **Combined PID histograms** (TPC + TOF):
  - **ptSpectrumPionCombinedPID**
  - **ptSpectrumKaonCombinedPID**
  - **ptSpectrumProtonCombinedPID**
  - **ptSpectrumElectronCombinedPID**
- Particle-specific histograms:
  - **ptHistogramPion/Kaon/Proton/Electron**: Reconstructed $p_T$
  - **dEdxPion/Kaon/Proton/Electron**: $dE/dx$ vs $p_T$ for identified particles

## Simulated (MC Truth) Histograms:
These histograms use **MC particle-level** information:

- **ptGeneratedPion/Kaon/Proton/Electron**: Generated $p_T$ of true MC particles
- **numberOfRecoCollisions**: Number of reconstructed collisions per MC collision
- **multiplicityCorrelation**: Track multiplicity correlation between reconstructed collisions

## Analysis Task Features:
The associated analysis task ([myExampleTaskPID.cxx](myExampleTaskPID.cxx)) includes:
1. Track selection with DCA cut (< 0.2 cm)
2. PID validation using:
   - TPC signal ($dE/dx$) vs momentum
   - $n\sigma$ separation for pions/kaons/protons/electrons
3. MC truth matching for:
   - $p_T$ resolution studies
   - Physical primary identification ($|y| < 0.5$)
4. Collision multiplicity studies

## Summary of Classification:
| Reconstruction Histograms          | MC Truth Histograms         |
|------------------------------------|-----------------------------|
| `eventCounter`                     | `ptGenerated*`              |
| `etaHistogram`                     | `numberOfRecoCollisions`    |
| `ptHistogram*`                     | `multiplicityCorrelation`   |
| `dEdx*`                            |                             |
| `tpcDedxVsPt`                      |                             |
| `nSigma*VsPt`                      |                             |

## Key Notes:
- The analysis task includes advanced features such as PID validation using both $dE/dx$ and $n\sigma$-based methods.
- Histograms like `ptResolution` require MC truth matching to compare reconstructed and generated particle properties.
- The `multiplicityCorrelation` histogram is used to study relationships between reconstructed collisions in events with multiple interactions.
