# CRISPR-Cas9 in T Cells — Molecular Editing Dynamics Model

This project formulates CRISPR-Cas9 genome editing in T cells as a **molecular-scale dynamical system**, where sequence recognition, binding stability, chromatin accessibility, and stochastic repair jointly determine editing outcomes.

Rather than treating editing as a binary event, the pipeline models the full probabilistic landscape of cleavage and repair. Biophysical priors are combined with statistical learning to estimate calibrated on- and off-target cutting probabilities, predict indel distributions, and propagate uncertainty through end-to-end scenario simulations.

CRISPR is treated analogously to an enzyme–substrate system: guide–DNA interactions define an energetic landscape, reaction pathways shape outcome distributions, and population-level behavior emerges from molecular-scale kinetics.

By integrating sequence context, epigenetic state, and probabilistic simulation, the framework provides a unified computational model of genome engineering decisions with standardized outputs designed for downstream systems-level modeling.

The broader aim is to connect molecular biophysics and modern learning methods into a rigorous computational framework for high-resolution modeling of biological interventions.
