# FlexiBLEPlugin
A plugin to calculate FlexiBLE force
##Intallation Instruction
A little heads-up:
The plugin replies on Conda Environment to install. 
The key software package will be:
  OpenMM        8.0.0
  CMake         3.27.6
  gcc_linux-64 >11
CUDA and OpenCL are not supported for this plugin. 
Python APIs not configured yet. 

1. Clone this repository to local machine. 
2. Use CCMake to configure, the CMAKE_INSTALL_PREFIX and OPENMM_DIR should be the same as your CONDA_PREFIX.
3. Generate when you make sure the paths are right.
4. Type "make install". 
5. If you wanna test FlexiBLE with molecular dynamics, alter the parameters of those tests within "platform/references/tests". The sample initial coordination and initial velocity files are within the "test" directory. 

This plugin is based on Zhuofan Shen and William J. Glover's work published on J. Chem. Phy., the paper "https://doi.org/10.1063/5.0067855" must be cited if using this plugin in research. 



