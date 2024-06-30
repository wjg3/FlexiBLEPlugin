# FlexiBLEPlugin
An [OpenMM](https://github.com/openmm/openmm) plugin that implements the FlexiBLE QM/MM boundary potential [[10.1063/5.0067855](https://doi.org/10.1063/5.0067855)]

## Intallation Instruction
A little heads-up:  
The plugin relies on several packages to install, will be easy to configure with Conda.  
The key software packages include: 

    OpenMM         8.0.0 and higher (7.3.0 works when changing LangevinMiddleIntegrator to LangevinIntegrator in tests) 
    CMake          3.27.6  
    gcc_linux-64 >=11  

CUDA and OpenCL are not supported for this plugin.  
Python APIs not configured yet. 

1. Clone this repository to local machine. 
2. Use CCMake to configure, the CMAKE_INSTALL_PREFIX and OPENMM_DIR should be the same as your CONDA_PREFIX.
3. Generate when you make sure the paths are right.
4. Type "make install". 
5. If you want to test FlexiBLE with molecular dynamics, alter the parameters of the tests within "platform/references/tests". The sample initial coordination and initial velocity files are within the "test" directory. 

## Citation info
The following must be cited if using this plugin in published research: 

1. Z. Shen and W. J. Glover, "Flexible boundary layer using exchange for embedding theories. I. Theory and implementation" *J. Chem. Phys.* **155**, 224112 (2021)
2. FlexibLEPlugin. An OpenMM plugin implementing the FlexiBLE QM/MM boundary potential. https://github.com/KaiChen-lr/FlexiBLEPlugin



