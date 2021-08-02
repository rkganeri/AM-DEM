# AM-DEM

A discrete element method (DEM) code for simulating powder deposition during additive manufacturing (AM).
A full description of the governing equations and basic methodology is provided in the following paper:

Ganeriwala, R., & Zohdi, T. I. (2016). A coupled discrete element-finite difference model of selective 
laser sintering. Granular Matter, 18(2). https://doi.org/10.1007/s10035-016-0626-0

A prior, open-access paper describing the basics (without the modified Hertzian contact model) is provided 
here:
Ganeriwala, R., & Zohdi, T. I. (2014). Multiphysics modeling and simulation of selective laser sintering 
manufacturing processes. In Procedia CIRP (Vol. 14). https://doi.org/10.1016/j.procir.2014.03.015

This code is based off some sloppy Fortran code I originally wrote in grad school to simulate the process.
The old code is contained in the original\_code folder, which was used to produce the results of the
2016 paper.

The new code within src was re-written in C++ and is capable of running on GPUs via the Kokkos offloading
library.  If Kokkos and this code are compiled with either CUDA or HIP the calculations can be performed
on either Nvidia or AMD GPUs, respectively.  However, if only CPUs are available when compiling, 
Kokkos still provides parallelism via OpenMP threading.  Note that my 2011 Macbook Pro laptop I've been 
running/testing these simulations on does not have any compute capable GPUs so all testing has been CPU-only. 
While I have coded this in a manner such that it should run correctly on GPUs (without unified memory)
as well, I have not tested that so there's a chance some small bugs are present.

Note that Kokkos automatically optimizes memory usage for multi-dimensional arrays (views) based upon the 
compilation configuration (e.g. LayoutRight memory for CPU-only or LayoutLeft for GPUs).  Additionally, using
Kokkos we are able to use either OpenMP or CUDA/HIP (depending on compilation setting) without the need for
any #ifdef's within the code!



Basic Compilation Instructions:

1) This repo utilizes submodules so when cloning the repo be sure to use the command:
```
git clone --recurse-submodules
```

If you have already cloned the repo without the submodules, you can add them later via: 
``` 
git submodule init
git submodule update
```


2) Build the third party libraries using the scripts provided in TPLs/TPL\_build\_scripts.

  a) kokkos

     Kokkos is used for offliading to GPUs or generating OpenMP code for execution on CPUs. To build Kokkos
     using CMake, enter the subdirectory TPLs/kokkos and create a build_* subdirectory therein, where
     the subdirectory should be named to indicate the build configuration (see names of build scripts for
     example). Copy, the appropriate Kokkos build script from TPLs/TPL\_build\_scripts depending on the 
     desired configuration, or generate your own build command based upon your architecture and desired
     settings. A few templates scripts are provided but for more build options visit the Kokkos github 
     page directly. Note that unlike all following TPLs, the AM-DEM build finds the appropriate Kokkos
     build directory in the host-config file, as it is anticipated many different version of Kokkos may be
     pre-compiled and linked to depending on AM-DEM build type settings.

  b) fmt

     Fmt is a library used for ease of outputting text, especially for debug purposes. It contains functions
     for generating strings in a manner similar to the python ".format()" method. To compile fmt create
     a build subdirectory within TPLs/fmt. Then copy the build script from TPLs/TPL\_build\_scripts and execute
     it within the build directory. Note that fmt may already be built on some systems, and there is a header
     only version of the library which may also be used and won't require separate compilation.

  c) lean-vtk

     Lean-vtk is used for ease of generating VTK plot files. To build it simply create a build directory
     within TPLs/lean-vtk, copy the appropriate shell script from the TPL\_build\_scripts folder here, and
     then execute the script.

  d) CLI11
     
     CLI11 is a header only library used to parse the command line. No compilation is necessary, it is 
     automatically linked to in the CMakeLists.txt file in the home folder of the AM-DEM repo.

  e) blt

     BLT is a library of CMake wrapper functions/macros that is useful for building, linking, and testing via
     CMake. It is used to automatically link in GTest for unit testing and adds some simple "smoke" tests
     to verify basic functionality of e.g. GTest, OpenMP, Cuda, and/or MPI. Again, nothing needs to be done 
     for this library. Similar to fmt, lean-vtk, and CLI11 it will be found and configured within the main 
     CMakeLists.txt file.


3) Build the current code using CMake. Create a build directory within the home directory for this repo. Run
   CMake using the appropriate host-config file for your build type, desired settings, and compiler.
   A few template host-config files using the Clang 11.1 compiler on a Mac operating the High Sierra OS 
   (10.13.6) are provided. For other example host-config files look at the BLT github repo and corresponding
   documentation, especially if you are linking to CUDA via nvcc, as that can sometimes be a bit tricky.
   To actually compile the executable once all TPLs are built, run from the home directory of the repo e.g.:
   ```
   mkdir build_clang11_debug
   cd build_clang11_debug
   cmake -C ../host-configs/macos_10_13_x86_64_clang11_debug.cmake ..
   make -j
   ```

   To execute the test problems, which should all pass, simply run:
   ```
   make test
   ```

   

Basic Running Instructions:

Many global parameters are set in the file src/global\_settings.hpp.  The easiest place to modify size of 
bounding volume and particle size is within that file.  The only required argument for running the code is 
the number of particles via the -n flag.  600 particles were used in the 2016 paper.  Plot files are writtin 
in VTK format, and the results can be visualized using an open source vis package such as Visit or Paraview.



Feel free to shoot me a message with any questions, bug reports, or suggestions for improvement.

