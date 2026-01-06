## Dynamics of a bricklayer model: multi-walker realizations of true self-avoiding motion

Code for the paper [https://www.arxiv.org/abs/2510.00975](https://www.arxiv.org/abs/2510.00975)
# Introduction to `bricks.cc`

This C++ program simulates the **Bricklayer Model**, a multi-walker generalization of the true self-avoiding walk described in the accompanying paper. It models the collective dynamics of $N$ agents ("bricklayers") building a wall by moving based on local height gradients.

## Core Concepts

The simulation evolves the system in continuous time using a **Rejection-Free Kinetic Monte Carlo** (Gillespie-like) algorithm, efficiently managing event times using a sorted data structure.

### Key Data Structures

1.  **Event Queue (`multiset<shared_ptr<siteEvent>> events`)**:
    *   This is the heart of the simulation. It stores potential future jump events for all sites, sorted by time.
    *   It uses a C++ `std::multiset` (typically implemented as a Red-Black Tree), allowing for efficient insertion, deletion, and retrieval of the next event in $O(\log N)$ time.
    *   **Why `multiset`?** Unlike a standard priority queue, a multiset allows us to efficiently remove specific arbitrary events (when a neighbor's state changes, their future event time becomes invalid and must be updated).

2.  **Lattice State**:
    *   `vector<int> z`: Represents the local gradient of the wall height ($z_i = h_{i-1} - h_i$).
    *   `vector<int> rho`: Stores the number of bricklayers (density) at each site $i$.
    *   `vector<direction> dir`: Stores the direction (left/right) of the next scheduled event for each site.

3.  **Cross-Linking (`vector<shared_ptr<siteEvent>> wall`)**:
    *   A vector that maps site indices to their corresponding event objects in the `multiset`.
    *   This allows $O(1)$ access to a site's event data to quickly update it when its neighbors change.

## Algorithm Overview

The simulation proceeds by processing events one by one:

1.  **Find Next Event**: Extract the event with the smallest time $t_{min}$ from `events`.
2.  **Move Bricklayer**:
    *   The bricklayer at site $l$ moves left or right based on the pre-calculated decision.
    *   Update the physical state: decrement `rho[l]`, increment `rho[target]`, and update gradients `z`.
3.  **Update Rates**:
    *   The move changes the local configuration (gradients $z$) not just at $l$, but also for its neighbors.
    *   **Recalculate**: Call `reset_site()` for $l$ and its affected neighbors. This removes their old events from the `events` queue and inserts new ones based on the new local rates $r(z) \sim e^{\beta z}$.
4.  **Advance Time**: Global time is updated to $t_{min}$.

## Code Structure

*   `main()`: Handles command-line arguments and triggers either a single visualization run (`single_run`) or a data-gathering ensemble run (`multi_run`).
*   `reset_site(j, time)`: The "engine" of the local updates. computes the Poisson rates for site $j$, determines the next event time/direction, and updates the `events` queue.
*   `manage_event(l, time)`: Executes the physical move (updating `z`, `rho`) and triggers `reset_site` for affected sites.
*   `rate(z)`: The exponential rate function driving the local gradient descent dynamics.

## Output Files

The program generates several data files which correspond to the results discussed in `nature.tex`:

*   **`width.dat`**: Evolution of the width (standard deviation) of the wall height distribution over time. Used to verify the $t^{2/3}$ scaling.
*   **`multi.dat`**: The averaged density profile $\rho(x)$ and height profile $h(x)$ from ensemble runs (`multi_run`). This data produces the parabolic/triangular scaling plots.
*   **`fourier.dat`**: Fourier components of the density, used to analyze traveling waves and oscillations.
*   **`final.dat`**: Final density distribution of builders.
*   **`zcor.dat`**: Correlations of the height gradient.


### Command Line Arguments

The program accepts up to 5 arguments:
1. `peak`: Number of builders at center if `dopeak` is true.
2. `nrun`: Multiplier for number of runs (nrun = 1000 * arg).
3. `tstep`: Time step size.
4. `N`: System size.
5. `ns`: Number of steps in single run.

Example:
```
./bricks 300 80 10 1024 1024


### Parameters (in namespace `Param`)

- `N`: System size (default 1024).
- `beta`: Parameter for rate function (default 0.4).
- `sweep_length`: Number of events per sweep (default N).
- `nrun`: Number of runs for multi-run mode (default 80000).
- `ns`: Number of steps in single run (default 1024).
- `tstep`: Time step size (default 10).
- `dopeak`: If true, places many builders at the center; else random placement.
- `GRAPH`: If true, enables animation.
