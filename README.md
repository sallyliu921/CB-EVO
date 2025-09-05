# CB-EVO
Logic synthesis sequence generation project

## Build the project
To get started with CB-EVO, clone the repository and build the project using the following commands:
```bash
git clone https://github.com/sallyliu921/CB-EVO
./build.sh
```

## Prerequisites

Ensure the following dependencies are installed and properly configured before building CB-EVO:

- **C++17 Compiler:**
  - Tested with GCC 9.4.0

- **OpenMP:**
  - Version tested: 4.0.1
  - Download and install commands:
    ```bash
    wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz
    tar -xvf openmpi-4.0.1.tar.gz
    cd openmpi-4.0.1
    ./configure
    make -j12 all install
    ```

- **Eigen3:**
  - Ensure Eigen3 is installed in your system. It can typically be installed via your package manager.

- **Yosys:**
  - Installation method:
    ```bash
    # Download OSS CAD Suite from GitHub
    wget <URL_to_OSS_CAD_Suite_release>
    tar -xzf <OSS_CAD_Suite_archive>
    # Add to your PATH
    echo 'export PATH="<extracted_location>/oss-cad-suite/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    ```
  - Replace `<URL_to_OSS_CAD_Suite_release>` and `<OSS_CAD_Suite_archive>` with the actual URL and downloaded archive name from [OSS CAD Suite Build](https://github.com/YosysHQ/oss-cad-suite-build/) appropriate for your system.
  
## Usage
Navigate to the build directory and run the program with specified command line arguments:
```bash
cd build/
./cb-evo -i <filename> -t <target> -s <sequence length> -e <exp_num> [-f]

# Suggested Example
./cb-evo -i ../benchmarks/EPFL/max.blif -t 2 -s 25
./cb-evo -i ../benchmarks/EPFL/max.blif -t 2
./cb-evo -i ../benchmarks/VTR8.0/bfly.blif -t 0 -s 42 -e 1 

```

### Command Line Arguments
- `-i <filename>`: **Required**. Specifies the input file.
- `-t <target>`: **Required**. Sets the optimization target for the synthesis sequence. Available options:
  - `0`: Optimize for AIG nodes (minimizing the number of nodes in And-Inverter Graphs)
  - `1`: Optimize for AIG nodes * level (minimizing both node count and circuit depth)
  - `2`: Optimize for LUTs (minimizing the number of Look-Up Tables)
  - `3`: Optimize for LUTs * level (minimizing both LUT count and circuit depth)
- `-s <sequence length>`: Sets the total sequence length for optimization. Range: 1-400. 
  - If not specified (default=0): 
    - Bandit algorithm will be disabled
    - Evo algorithm will dynamically extend the sequence based on case size, using runtime as the termination condition.
- `-e <exp_num>`: Sets the experiment mode. Range: 0-1.
  - `0` (default): Normal mode
    - For target 2 or 3 (LUT optimization), uses FPGA mapping command 'if -a -K 6'
  - `1`: FlowTune comparison mode
    - For target 0 or 1 (AIG optimization): Adds "ifraig; dch -f;" every 14 commands
    - For target 2 or 3 (LUT optimization): Adds "ifraig;scorr;dc2;strash;dch -f;if -K 6;mfs2;lutpack -S 1;" every 14 commands
- `-f`: Enable ccirc and yosys for additional static feature extraction in bandit model (default=0).

### Reference
```shell
@inproceedings{liu2024cbtune,
  title={CBTune: Contextual Bandit Tuning for Logic Synthesis},
  author={Liu, Fangzhou and Pei, Zehua and Yu, Ziyang and Zheng, Haisheng and He, Zhuolun and Chen, Tinghuan and Yu, Bei},
  booktitle={2024 Design, Automation \& Test in Europe Conference \& Exhibition (DATE)},
  pages={1--6},
  year={2024},
  organization={IEEE}
}
```