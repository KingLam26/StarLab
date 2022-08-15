# StarLab Optics

    .
    ├── beamClass.py            # support module, allows for the definition and manipulation of light beams
    ├── functions.py            # support module, handles line and circle geometry
    ├── functions.py            # support module, handles line and circle geometry
    ├── functions.py            # support module, handles line and circle geometry
    ├── functions.py            # support module, handles line and circle geometry

├── StarLab/
  ├── 
  ├── 
  ├── main.py           LSDS
  ├── main4.py          Horiba
  ├── main5.py          EQ-77
  ├── README.md

While it is certain that one could find existing optics simulation software that are more powerful and comprehensive than what is hastily put together here, they are often too complicated for our purpose, and the effort needed to properly understand them to obtain meaningful results is non-trivial. 

StarLab Optics seeks to use simple geometry and basic principles of refraction and reflection to quickly generate a baseline model that can predict how a finite beam of light behaves when manipulated using simple optics like lens, windows and mirrors etc, in an ideal world. It is always possible to slap on increasingly complex models to provide predictions with greater accuracy and precision. Nevertheless, it is important to bear in mind the practicalities of setting up an optics experiment, and to leave ample room for fine-tuning onsite, rather than design a fixed experimental setup.

An optimization module from Scipy is used to solve for the positions of certain pieces of optics while targetting a particular desired parameter.

For definitions of the syntax in the mainX.py files used to describe the beam positions and angles, please kindly refer to Star Lab Update 2.2, dated 15 Aug 2022.
