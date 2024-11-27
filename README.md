# StarLab Optics

The objective is to design an experimental stellar simulator called StarLab to study prebiotic chemical reactions under various planetary atmospheres. One key design requirement is simultaneous broadband illumination (UV, VIS, IR) focused onto either a 10-mm-sided quartz cuvette sample holder or 50-mm-diameter watch glass, with provisions to minimise photon loss from source to sample.

When it comes to optics design, one can certainly find existing optics simulation software that are more powerful and comprehensive than what is hastily put together here. But the latter are often too complicated (and expensive) for our purpose, and the effort needed to properly understand them to obtain meaningful results is non-trivial.

StarLab Optics seeks to use simple geometry and basic principles of refraction and reflection to quickly generate a baseline model that can predict how a finite beam of light behaves when manipulated using simple optics like lens, windows and mirrors etc, in an ideal world. It is always possible to slap on increasingly complex models to provide predictions with greater accuracy and precision. Nevertheless, it is important to bear in mind the practicalities of setting up an optics experiment, and to leave ample room for fine-tuning onsite, rather than design a fixed experimental setup.

An optimization module from Scipy is used to solve for the positions of certain pieces of optics while targetting a particular desired parameter.

# Repository

```
    .
    ├── StarLab/
        ├── beamClass.py            # support module, allows for the definition and manipulation of light beams
        ├── functions.py            # support module, handles line and circle geometry
        ├── main-1.py               # LSDS
        ├── main-4.py               # Horiba
        ├── main-5.py               # EQ-77
```
