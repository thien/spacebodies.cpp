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

## Consistency & Convergence

<!-- 
Study two particles. Choose their velocity and initial position such
that they directly collide with each other. Create a table where you document the
resulting/merged particles position. Study this position for various time step size choices
(switching off the adaptive time step size choice from Step 1 question 3 above) and compare
these positions to runs where you switch on the feature from Step 1 question 3. Can you
uncover the convergence order of your scheme? Compare the obtained accuracy also to cost
in terms of time step count.
-->

The table below describes the information related to the two bodies used for the collision experiment.

| Body # |  $s(x)$ |  $s(y)$ |  $s(z)$ |  $v(x)$ |  $v(y)$ |  $v(z)$ | Mass |
|--------+---------+---------+---------+---------+---------+---------+------|
| 1 | +0.1 | +0.1 | +0.1 | -2.0 | -2.0 | -2.0 | $1e^{-11}$ |
| 2 | -1.0 | -1.0 | -1.0 | +2.0 | +2.0 | +2.0 | $1e^{-11}$ |

The bodies collide at $(x,y,z) = (0.155,0.155,0.155)$ at time $t=0.21$. The error is calculated between the position of the error of the new item. The following table shows the timestep used to calculate the collision, and the error value left over. We calculate the error of only one dimension; $x$, as the others would follow a very similar pattern.

| Timestep        | Error at $x$     | Collision? |
|-----------------+------------------+------------|
|0.001 (Adaptive) | -0.00000000351134| yes        |
|0.000000005      | -0.00000000353880| yes        |
|0.0000000025     | -0.00000000337228| yes        |
|0.00000000125    | -0.00000000481655| yes        |
|0.000000000625   | -0.00000000509412| yes        |
|0.0000000003125  | -0.00000000540662| yes        |
|0.00000000015625 | -0.00000000649210| yes        |

Our convergence scheme 

The adaptive timestep works suitably and utilises less timesteps than $ts=5e^-9$ by a significant margin, whilst having relatively comporable error residues. <!-- Will need to talk about calculating at earlier timesteps to see whether they collide.-->

## Complexity
<!-- 
Run your code for 10, 100, 1,000, 10,000 particles placed randomly in space.
Derive the runtime complexity of the code and compare it to your experimental data.
-->

Under the assumption that the timestep and time limit is fixed, the most dominant function `updateBodies()` which utilises a nested loop that iterates through the number of bodies intiated. For each iteration, a force for a given body is calculated by comparing its position against every other body in space. This results in `updateBodies()` to run in $O(n^2)$. 

Procedures have been taken to reduce the constant; Each body only calculates its force against bodies that precede them in the order of initiation i.e. Body $2$ calculates force from Body $1$ and Body $0$ whereas Body $3$ calculates from $0,1$ and $2$. A pseudocode describes this method:

    for (i=0; i<n):
      for (j=0; j<i):
        if (i != j):
          distance = distance between body i and body j

          m = i.mass*j.mass/distance/distance/distance;

          for (k = 0; k < 3; k++):
            i.force[k] += (i.x[k]-j.x[k]) * m;
            j.force[k] += (i.x[k]-j.x[k]) * m;

Whilst `updateBodies()` would continue to run in $O(n^2)$, the hidden constant would be drastically reduced to a factor of $\frac{1}{2}$ of the original number of calculations needed.

## Statistics
<!-- 
Extend your code such that it keeps track of the number of bodies over time.
Create a plot that shows how the total number of particles decreases over simulation time as
particles merge.
-->

# Scaling Experiments

<!-- Repeat the experiments from Step 2 to ensure that your modifications did not break the code. From
hereon, create scaling plots for 10-10,000 particles. You are strongly encouraged to use a University
machine for your plots that has at least 4 cores, i.e. you present a scaling plot than spans at least
1,2,3 and 4 cores. If you have a more powerful machine at home, you are free to use this machine.
Clarify explicitly in your report the machine specifica. -->

The machine used consists of Durham's MIRA machine, which is a 128 core intel xeon distributed system. At runtime, the program consumes less than a megabyte of memory. 

# Questions

1. How does the scalability for very brief simulation runs depend on the total particle count?

2. Calibrate Gustafson’s law to your setup and discuss the outcome. Take your considerations on
the algorithm complexity into account.

3. How does the parallel efficiency change over time if you study a long-running simulation?

# Distributed Memory Simulation
<!-- 
Design a strategy how to parallelise your code with MPI. No implementation is required, i.e. it is a
gedankenexperiment. -->
