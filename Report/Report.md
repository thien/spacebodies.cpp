---
title: "Space Bodies Assignment"
author: cmkv68
date: February 19, 2018
geometry: "left=2cm,right=2cm,top=2cm,bottom=2cm"
output: pdf_document
---

<!--
The report may not exceed two pages (incl. pictures/figures) for the Steps 2,4 and 5, 
i.e. it may not exceed six pages in total.
-->

# Numerical Experiments

**Consistency/convergence**: Study two particles. Choose their velocity and initial position such
that they directly collide with each other. Create a table where you document the
resulting/merged particles position. Study this position for various time step size choices
(switching off the adaptive time step size choice from Step 1 question 3 above) and compare
these positions to runs where you switch on the feature from Step 1 question 3. Can you
uncover the convergence order of your scheme? Compare the obtained accuracy also to cost
in terms of time step count.

**Complexity**: Run your code for 10, 100, 1,000, 10,000 particles placed randomly in space.
Derive the runtime complexity of the code and compare it to your experimental data.

The complexity of the code runs in $$O(n^2)$$.

**Statistics**: Extend your code such that it keeps track of the number of bodies over time.
Create a plot that shows how the total number of particles decreases over simulation time as
particles merge.

# Scaling Experiments

Repeat the experiments from Step 2 to ensure that your modifications did not break the code. From
hereon, create scaling plots for 10-10,000 particles. You are strongly encouraged to use a University
machine for your plots that has at least 4 cores, i.e. you present a scaling plot than spans at least
1,2,3 and 4 cores. If you have a more powerful machine at home, you are free to use this machine.
Clarify explicitly in your report the machine specifica.

# Questions

1. How does the scalability for very brief simulation runs depend on the total particle count?
2. Calibrate Gustafson’s law to your setup and discuss the outcome. Take your considerations on
the algorithm complexity into account.
3. How does the parallel efficiency change over time if you study a long-running simulation?

# Distributed Memory Simulation

Design a strategy how to parallelise your code with MPI. No implementation is required, i.e. it is a
gedankenexperiment.
