# CHC: Mining and Visualizing Maximum Cohesive Hyper Clique

This repository provides a C++ implementation for the paper:

**Mining and Visualizing Maximum Cohesive Hyper Clique**  
Xiaoyu Leng, Hongchao Qin, Rong-Hua Li (Beijing Institute of Technology)

The codebase focuses on mining **maximum cohesive hyper cliques (CHC)** in hypergraphs, and includes multiple implementations/indexing strategies, including basic construction, heap-accelerated construction, approximate indexing, and non-index baselines.

## 1. Problem Definition (CHC)

Given a hypergraph \(H=(V,E)\), a vertex set \(C \subseteq V\) is an **η-cohesive hyper clique (η-CHC)** if **every pair of distinct vertices** in \(C\) co-occurs in at least **η common hyperedges**.  
This generalizes classical graph cliques: when all hyperedges have size 2, an η-CHC with η=1 reduces to a graph clique.

The goal of this artifact is to compute **maximum** CHCs (largest size) efficiently, using indexing and pruning strategies described in the paper.

## 2. Repository Structure

```
.
├── src/
│   ├── max_basic.cpp
│   ├── max_heap.cpp
│   ├── max_noindex.cpp
│   ├── max_approx.cpp
│   ├── max_heap_eta.cpp
│   └── max_approx_eta.cpp
├── dataset/
│   └── (hypergraph datasets or download instructions)
└── README.md
```


## 3. Requirements

- macOS / Linux
- A C++17 compiler:
  - `g++ >= 9` or `clang++ >= 10` (C++17 supported)

No external dependencies are expected beyond the C++ standard library unless your code explicitly requires them (please update this section if you use third-party libs).

## 4. Build

### Option A: Build with `g++` / `clang++` (simple)

From repository root:

```bash
mkdir -p bin

# Build each variant (adjust file paths if you keep sources in root)
c++ -O3 -std=c++17 -o bin/max_basic        src/max_basic.cpp
c++ -O3 -std=c++17 -o bin/max_heap         src/max_heap.cpp
c++ -O3 -std=c++17 -o bin/max_noindex      src/max_noindex.cpp
c++ -O3 -std=c++17 -o bin/max_approx       src/max_approx.cpp
c++ -O3 -std=c++17 -o bin/max_heap_eta     src/max_heap_eta.cpp
c++ -O3 -std=c++17 -o bin/max_approx_eta   src/max_approx_eta.cpp
```



## 5. Input Data Format

Expected dataset format (common hypergraph incidence format)

Each hyperedge is represented by a list of vertex IDs (integers) on one line:
v1 v2 v3 ...

- One hyperedge per line.
- Vertex IDs should be non-negative integers.
- Repeated vertex IDs within a line should not appear (treat as a set).

## 6. Usage

This repository includes multiple executables, corresponding to different strategies described in the paper:

- `max_basic`: basic CHC-index construction and maximum CHC
- `max_heap`: heap-accelerated CHC-index construction
- `max_noindex`: baseline without index (for comparison)
- `max_approx`: approximate index for faster runtime with comparable solution quality
- `*_eta`: variants emphasizing η-related optimizations / η-peak interval logic (as in the paper)



