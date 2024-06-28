Main repository for the reinforcement-learning based path-planning Summer Internship Project supervised by Dr. Zhongguo Li at the University of Manchester.

Present aims of the project seeks to convert MATLAB example code into Python to leverage usage of PyTorch.

The code has been fully translated through to Python, but still requires debugging to remove errors and check correct operation.

Attached below are screenshots of replicated operation thusfar:

__Particle Filter initial (x,y) scatterplot & 3D Gaussian Plume__

![image](https://github.com/b-kirk/Reinforcement-Learning-based-Path-Planning/blob/main/Screenshots/Particle_Filter_3D.png)

__Histogram of generated Gamma Distribution__

![image](https://github.com/b-kirk/Reinforcement-Learning-based-Path-Planning/blob/main/Screenshots/Gamma_Dist.png)

_Graphs generated using the Plotly Python library_

The steps of the project following complete translation are to integrate other concentration curves than the Gaussian Plume (e.g. generating a point-source radioactive emitter), as well as general improvements to path-length and other parameters of the algorithm.
