# TrackMPD

TrackMPD is a three-dimensional particle-tracking model for the transport of marine plastic debris in oceans and coastal systems. The power of TrackMPD lies in: (1) its compatibility with diverse formats of current-velocity inputs; and (2) its ability to extend the Lagrangian modelling of advection-diffusion by adding more-complex and realistic particle behaviours and physical processes, which can either be included or excluded depending on the application. At present, TrackMPD can include beaching, washing-off, degradation, biofouling, sinking, and deposition. In particular, sinking and deposition depend on particle behaviour, which relies on the particle density, size, shape, fouling state, and degradation state. The model can incorporate new processes and behaviours, and change the implementation of already existing ones, with new experimental findings or particular applications.

TrackMPD has thus a structured and coherent modelling framework to satisfy the criteria of flexibility, extendability, and interchangeability. This framework is based on the Particle Tracking and Analysis Toolbox (PaTATO, Fredj et al., 2016) for the advection and diffusion processes. It can use velocity data from various sources, such as different ocean general circulation models (OGCM: e.g. POM, FVCOM) and satellite observations, and can compute forward and backward trajectories in two or three dimensions. The model consists of a set of coupled and mutually interacting modules. Modules are independent functions or classes that define behaviours, read the inputs from a certain source, implement a given physical process, or perform auxiliary tasks such as creating outputs. This allows independent development of modules that can be easily added to the model without the need to change the other modules. 

From Jalón-Rojas et al. (2019). Refer to this publication for further details.


# Read more and getting started

Several options are provided to learn to work with TrackMPD.

- A getting-started tutorial — see folder Manuals.
- Example input files linked to (scientific) publications or other applications of the model — see folder Examples:
	* POM Application: Jalón-Rojas, I., Wang, X.-H., and Fredj, E.: Technical note: On the importance of a three-dimensional approach for modelling the transport of neustic microplastics, Ocean Sci., 15, 717-724, https://doi.org/10.5194/os-15-717-2019, 2019.
	* FVCOM Application: developed by Zhixin Cheng. 
- Extensive manuals for TrackMPD will be released with version 2.

Our paper “A 3D numerical model to Track Marine Plastic Debris (TrackMPD): Sensitivity of microplastic trajectories and fates to particle dynamical properties and physical processes” is published and available via [this link](https://www.sciencedirect.com/science/article/pii/S0025326X19301523).

TrackMPD v.2 will be release soon, together with a new publication.


# Programming language and prerequisites

TrackMPD is developed for MATLAB. We recommend installing the [M_MAP toolbox](https://www.eoas.ubc.ca/~rich/map.html#9._Zoom_in_on_Prince_Edward_Island_to_co) or the function “m_fdist” of this toolbox before using TrackMPD.
 
Running TrackMPD requires very little knowledge of MATLAB or programming in general. Users might, however, need to use TrackMPD with a specific OGCM, which requires basic knowledge of MATLAB. Authors are happy to help users to make TrackMPD compatible with different OGCMs.


# License and terms of use

When using TrackMPD in any scientific publication, technical report or otherwise formal writing, please cite our paper:
Jalón-Rojas, I., Wang, X.H., Fredj, E., 2019. “A 3D numerical model to Track Marine Plastic Debris (TrackMPD): Sensitivity of microplastic trajectories and fates to particle dynamical properties and physical processes”. Marine Pollution Bulletin, 141, 256-272.

The TrackMPD code is licensed under GPL (GNU General Public License). In summary, this means that the code is open source and may be used freely for non-commercial and commercial purposes. Any alterations to the TrackMPD source code or new modules must be licensed under GPL as well. See the LICENSE.md.


# Authors

* Erick Fredj: Main developer, responsible for the TrackMPD v.2 core, and the advection and dispersion modules (erick.fredj@gmail.com).
* Isabel Jalón-Rojas: Main developer, responsible for the TrackMPD v.1 core, and the behaviour, settling and washing-off modules (ijalonrojas@gmail.com).
* Xiao-Hua Wang: Contributor. 


# Contributing

Contributions to TrackMPD are welcome. Contact authors via email for more information.
