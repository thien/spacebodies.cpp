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

The bodies collide at $(x,y,z) = (0.155,0.155,0.155)$ at time $t=0.275$. The error is calculated through calculating the difference between the position of one body and the other body upon collision. The following table shows the timestep used to calculate the collision, and the error value left over. We calculate the error of only one dimension; $x$, as the other dimensions ($y,z$) would follow the same values.


| Timestep        | Error    | Position of $x$ |
|-----------------+------------------+------------|
|0.0001 (Adaptive) | 0.000199998| -0.4498        |
|0.000078125000 | 0.000156250000 | -0.449844 |
|0.000039062500 | 0.000078125000 | -0.449922 |
|0.000019531300 | 0.000039062500 | -0.449961 |
|0.000010000000 | 0.000020000000 | -0.44998 |
|0.000003906250 | 0.000007812500 | -0.449992 |
|0.000001953130 | 0.000003906250 | -0.449996 |
|0.000001000000 | 0.000001999990 | -0.449998 |
|0.000000610352 | 0.000001220720 | -0.449999 |
|0.000000305176 | 0.000000610365 | -0.449999 |
|0.000000152588 | 0.000000305098 | -0.45 |


![A chart showing the the timestep used against the error produced. The chart shows a clear convergence.](chart.png)


<!-- 
COLLISION @ t=0.275, Timestep: 0.001: Error: 0.002, Location: -0.448, Ratio: 0.5
COLLISION @ t=0.275, Timestep: 0.0001: Error: 0.0002, Location: -0.4498, Ratio: 0.5
COLLISION @ t=0.275, Timestep: 1e-05: Error: 2e-05, Location: -0.44998, Ratio: 0.5
COLLISION @ t=0.275, Timestep: 1e-06: Error: 1.99999e-06, Location: -0.449998, Ratio: 0.500003
COLLISION @ t=0.275, Timestep: 1e-07: Error: 2.00017e-07, Location: -0.45, Ratio: 0.499958
COLLISION @ t=0.275, Timestep: 1e-08: Error: 1.97119e-08, Location: -0.45, Ratio: 0.507308
COLLISION @ t=0.275, Timestep: 1e-09: Error: 7.58933e-09, Location: -0.45, Ratio: 0.131764
COLLISION @ t=0.275, Timestep: 1e-10: Error: 1.51078e-08, Location: -0.45, Ratio: 0.0066191 -->

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

  The overhead involved in initiating a large number of particles for a set of parallel processors may take more time than the simulation itself. 

2. Calibrate Gustafson’s law to your setup and discuss the outcome. Take your considerations on the algorithm complexity into account.

  Gusafson's Law: 
  Gustafson estimated the speedup S gained by using N processors (instead of just one) for a task with a serial fractions (which does not benefit from parallelism) as follows:
  S=N+(1-N)s

3. How does the parallel efficiency change over time if you study a long-running simulation?

  There

# Distributed Memory Simulation
<!-- 
Design a strategy how to parallelise your code with MPI. No implementation is required, i.e. it is a
gedankenexperiment. -->


<!-- Mark Scheme

- All three features of Step 1 are implemented correctly.
- All three questions from Step 2 are answered correctly.
- Quality of presentation of simulation outputs for Step 2 (graphs, snapshots and notably videos).
- The code is parallelised including the statistics on the total number of objects.
- All three questions from Step 4 are answered correctly.
- A clear concept for MPI parallelisation is presented.
- Quality of writing and quality of presentation of simulation outputs (graphs, snapshots).

 -->