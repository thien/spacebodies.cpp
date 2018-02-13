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
<!-- 40 marks -->
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

| Timestep $h$ |  Error $\overline{u}_{h}$|  $x_{a}$ |  $x_{b}$ |  Ratio |  Steps |  Range |
|------+------+------+------+-----+------+------+------|
|Adaptive|0.000001997990|-0.449998|-0.450002|0.500503|1716457|0.000003994500|
|$10^{-6}/2^1$|0.000001999990|-0.449998|-0.450002|0.500003|275000|0.000003998500|
|$10^{-6}/2^2$|0.000001000020|-0.449999|-0.450001|0.499992|550000|0.000001998520|
|$10^{-6}/2^3$|0.000000499961|-0.45|-0.45|0.500039|1100000|0.000000998475|
|$10^{-6}/2^4$|0.000000250068|-0.45|-0.45|0.499865|2200000|0.000000498564|
|$10^{-6}/2^5$|0.000000125078|-0.45|-0.45|0.499688|4400000|0.000000248608|
|$10^{-6}/2^6$|0.000000062134|-0.45|-0.45|0.502946|8800000|0.000000123183|
|$10^{-6}/2^7$|0.000000031729|-0.45|-0.45|0.492452|17600000|0.000000061368|
|$10^{-6}/2^8$|0.000000014415|-0.45|-0.45|0.541966|35200000|0.000000028684|
|$10^{-6}/2^9$|0.000000006426|-0.45|-0.45|0.607915|70400000|0.000000012311|
|$10^{-6}/2^{10}$|0.000000002519|-0.45|-0.45|0.775391|140800000|0.000000004178|


(0.000000014415 - 0.000000006426) / (0.000000006426 - 0.000000002519)

The convergence order is $0.7152957794....$

![A chart showing the the timestep used against the error produced. The chart shows a clear convergence.](chart.png)

<!-- The numerical method is of order $p$ the there exists a number $C$ indepdendent of $h$ such that $ |u˜h − u| ≤ Chp$; at least for a sufficently small $h$.  -->

---

The adaptive timestep works suitably and utilises less timesteps than $ts=5e^-9$ by a significant margin, whilst having relatively comporable error residues. <!-- Will need to talk about calculating at earlier timesteps to see whether they collide.-->

## Complexity
<!-- 
Run your code for 10, 100, 1,000, 10,000 particles placed randomly in space.
Derive the runtime complexity of the code and compare it to your experimental data.
-->

Under the assumption that the timestep and time limit is fixed, the most dominant function `updateBodies()` which utilises a nested loop that iterates through the number of bodies initiated. For each iteration, a force for a given body is calculated by comparing its position against every other body in space. This results in `updateBodies()` to run in $O(n^2)$. 

Procedures have been taken to reduce the constant; Each body only calculates its force against bodies that precede them in the order of initiation i.e. Body $2$ calculates force from Body $1$ and Body $0$ whereas Body $3$ calculates from $0,1$ and $2$. 

<!-- A pseudocode describes this method:

    for (i=0; i<n):
      for (j=0; j<i):
        if (i != j):
          distance = distance between body i and body j

          m = i.mass*j.mass/distance/distance/distance;

          for (k = 0; k < 3; k++):
            i.force[k] += (i.x[k]-j.x[k]) * m;
            j.force[k] += (i.x[k]-j.x[k]) * m; -->

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

The machine used consists of a `Intel i7 3770k` processor at a 3.7Ghz clock speed, powering 4 cores and 8 threads. It utilises 32GB of memory and the storage consists of a SSD hooked up via SATA3. It is using a fresh installation of Ubuntu and has no other major processes running.

# Questions
<!-- 30 marks -->

1. __How does the scalability for very brief simulation runs depend on the total particle count?__

  There is a lot of overhead involved in initiating a loop for a set of parallel processors, to the extent that _may take more time than the actual simulation itself._ This may include each processor initiating their own set of variables.

2. __Calibrate Gustafson’s law to your setup and discuss the outcome. Take your considerations on the algorithm complexity into account.__

  Gusafson's Law: 
  Gustafson estimated the speedup S gained by using N processors (instead of just one) for a task with a serial fractions (which does not benefit from parallelism) as follows:
  S=N+(1-N)s

    I It depends on whether you fix the problem size.
  I It hence depends on your purpose.
  I It is crucial to clarify assumptions a priori.
  I It is important to be aware of shortcomings.

3. __How does the parallel efficiency change over time if you study a long-running simulation?__

  - Parallelism on a single machine now depends on other factors
    - makes more heat
    - poor code
    - transistor noise
    - contingent on other factors on the machine

# Distributed Memory Simulation
<!-- 20 marks -->
<!-- 
Design a strategy how to parallelise your code with MPI. No implementation is required, i.e. it is a
gedankenexperiment. -->

## Assumptions and Setup

- MPI is SPMD (Single program; multiple data)
  - we check the ID to determine whether it is a master or a slave, we call master rank 0.
  - the master has its own set of functions, and so does the master

- Assume perfect zero latency theoretical model between ranks

- For sake of simplicity the master-slave model is adopted; the master rank does not perform any major computation; this is distributed to the slaves

- Every CPU in the MPI network will have a copy of the same code, but is allocated different sets of data in regards to the loop position they are to compute

## Operation

- master rank initiates a large array of bodies and `MPI_Bcast` to slaves
- whilst this utilises a lot of data transmission, it is necessary as when calculating the force for a given body, it compares its position against every other body in the space.

- for calculating the force
  - Each rank will perform segments of the loop via `MPI_Scatter` (contingent on the number of bodies; it may be better given a large enough set (talking 100000+) to propagate them manually as the MPI buffers could heavily stress)
    - they're passed the relevant information;
      - master rank uses `MPI_send()` to a slave; slave receives using `MPI_recv()`
      - verify data is received properly through checking the status of the received message (this applies to all messages transmitted; if it isn't a valid response then we should request the message again or have error handlers in place.)
      - slaves have their own local variables that they can use to compute
      - also `gets` relevant information from the global space such as the positions and velocities of bodies
    - once the slave rank has finished computing the result, being the force increment for a body, are sent back to the master via `MPI_Reduce`
  - will need to check that we have all the results!

  - they also calculate
    - non blocking communication between the slave and master to send/receive information about short distances between two bodies as we calculate the potential collision afterwards; we use an `immediate_send`.
  
- for collisions
  - contingent on the number of collisions occurred
  - this can be distributed to nodes but depends on how fast data can be transmitted to

- for adaptive timestepping
  - we can distribute the potential solution between nodes

- for updating the bodies
  - given enough bodies this can be distributed across the network too using a `MPI_Scatter` and `MPI_Reduce`

- Process is repeated accordingly

## Issues

<!-- Mark Scheme

- All three features of Step 1 are implemented correctly.
- All three questions from Step 2 are answered correctly.
- Quality of presentation of simulation outputs for Step 2 (graphs, snapshots and notably videos).
- The code is parallelised including the statistics on the total number of objects.
- All three questions from Step 4 are answered correctly.
- A clear concept for MPI parallelisation is presented.
- Quality of writing and quality of presentation of simulation outputs (graphs, snapshots).

 -->