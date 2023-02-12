## Overview

These are simulation code files related to our papers "Integrated Path Following and Collision Avoidance Using a Composite Vector Field" and "Guiding Vector Fields for Following Occluded Paths". You could use the following bibtex items to cite these papers.

```tex
@inproceedings{yao2019integrated,
	title={Integrated Path Following and Collision Avoidance Using a Composite Vector Field},
	author={Yao, Weijia and Lin, Bohuan and Cao, Ming},
	booktitle={2019 IEEE 58th Conference on Decision and Control (CDC)},
	pages={250--255},
	year={2019},
	organization={IEEE}
}

@ARTICLE{yao2022collision,
  author={Yao, Weijia and Lin, Bohuan and Anderson, Brian D. O. and Cao, Ming},
  journal={IEEE Transactions on Automatic Control}, 
  title={Guiding Vector Fields for Following Occluded Paths}, 
  year={2022},
  volume={67},
  number={8},
  pages={4091-4106},
  doi={10.1109/TAC.2022.3179215}
}
```



## Code usage

Please run *.slx or *\_ode.m files first to get simulation result, after which you can run *\_animate.m file to get animation of the data and *\_plot.m to plot the data.



## Files and descriptions

| File              | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| d_mixvf...        | path: circle; 5 obstacle: ellipses and a Cassini oval; revised bump functions |
| f_switch...       | path: ellipse; 1 obstacle: curved square; stable equilibrium; switching vector field; |
| i_3Dunicycle      | 3D; unicycle/single-integrator;                              |
| j_moving_unicycle | path: sin curve; 1 moving obs: ellipse; unicycle/single-integrator; |

​      

   





