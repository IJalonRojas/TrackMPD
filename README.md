# TrackMPD

TrackMPD is a three-dimensional particle-tracking model for the transport of marine plastic debris in oceans and coastal systems. The power of TrackMPD lies in: (1) its compatibility with diverse formats of current-velocity inputs; and (2) its ability to extend the Lagrangian modelling of advection-diffusion by adding more-complex and realistic particle behaviours and physical processes, which can either be included or excluded depending on the application. At present, TrackMPD can include beaching, washing-off, sinking, deposition, resuspension, and bed load. In particular, sinking, deposition and resuspension depend on particle behaviour, which relies on the particle density, size, shape, fouling state, and degradation state. The model can incorporate new processes and behaviours, and it is possible to easily modify the implementation of the already existing ones with new experimental findings or particular applications.

TrackMPD has thus a structured and coherent modelling framework to satisfy the criteria of flexibility, extendability, and interchangeability. It can use velocity data from various sources, such as different global, regional and coastal ocean models (e.g. POM, NEMO, FVCOM, TELEMAC, HYCOM, MARS) and satellite observations, and can compute forward and backward trajectories in two or three dimensions. The model consists of a set of coupled and mutually interacting modules. Modules are independent functions or classes that define behaviours, read the inputs from a certain source, implement a given physical process, or perform auxiliary tasks such as creating outputs. This allows the independent development of modules that can be easily added to the model without the need to change the other modules. 

Modified from Jalón-Rojas et al. (2019). Refer to this publication for further details.


# TrackMPD versions

- TrackMPD_v3:
     * Same reference system as hydrodynamic models
     * New parameterizations for resuspension
     * New beaching calculation 
     * Water level over the particle's position is now saved in the output file


- TrackMPD_v2:

     * Parallel computation of particle trajectories
     * Optimized advection and behaviour modules
     * New processes: deposition, resuspension, bedload
     * New parameterizations of settling velocities
     * Lower computational time

- TrackMPD_v1:

     * Processes: advection, dispersion, settling, beaching, washing-off
     * Dynamical behaviour as a fonction of the particle properties and fouled state.


# Read more and getting started

Several options are provided to learn to work with TrackMPD.

- Full tutorial for TrackMPD v2.3 and TrackMPD v3.0.

- A getting-started tutorial of TrackMPD v1.

- Example input files linked to (scientific) publications or other applications of the model — see folder Examples:
	
For TrackMPD v2 and v3:

* HYCOM Application in the Ría de Arousa (Spain)

* TELEMAC Application in the GAronne tidal River

- Our paper "A wave-resolving two-dimensional vertical Lagrangian approach to model microplastic transport in nearshore waters based on TrackMPD 3.0" is published and availbe [here.](https://gmd.copernicus.org/articles/18/319/2025/) 

For TrackMPD_v1:

* POM Application: Jalón-Rojas, I., Wang, X.-H., and Fredj, E. (2019). Technical note: On the importance of a three-dimensional approach for modelling the transport of neustic microplastics, Ocean Sci., 15, 717-724, https://doi.org/10.5194/os-15-717-2019.

* FVCOM Application: Cheng, Z., Jalón-Rojas, I., Wang, X.H., Liu, Y. (2020). Impacts of land reclamation on sediment transport and sedimentary environment in a macro-tidal estuary, Estuar. Coast. Shelf Sci., 221, 106861. doi: https://doi.org/10.1016/j.ecss.2020.106861.


- Our paper “A 3D numerical model to Track Marine Plastic Debris (TrackMPD): Sensitivity of microplastic trajectories and fates to particle dynamical properties and physical processes” is published and available via [this link](https://www.sciencedirect.com/science/article/pii/S0025326X19301523) or [this one] (http://isabeljalonrojas.com/wp-content/uploads/2019/09/2019_mpb_JalonRojasetal.pdf).



# Programming language and prerequisites

TrackMPD is developed for MATLAB. We recommend installing the [M_MAP toolbox](https://www.eoas.ubc.ca/~rich/map.html#9._Zoom_in_on_Prince_Edward_Island_to_co) or the function “m_fdist” of this toolbox before using TrackMPD.
 
TrackMPD_v2 and TrackmPD_v3 requires the installation of the Matlab Parallel Toolbox.

Running TrackMPD requires very little knowledge of MATLAB or programming in general. Users might, however, need to use TrackMPD with a specific hydrodynamic model, which requires more advanced knowledge of MATLAB. Authors are happy to help users to make TrackMPD compatible with different hydrodynamic models.


# License and terms of use

When using TrackMPD in any scientific publication, technical report or otherwise formal writing, please cite our papers:

Jalón-Rojas, I., Wang, X.H., Fredj, E., 2019. “A 3D numerical model to Track Marine Plastic Debris (TrackMPD): Sensitivity of microplastic trajectories and fates to particle dynamical properties and physical processes”. Marine Pollution Bulletin, 141, 256-272.



Jalón-Rojas, I., Sous, D., and Marieu, V., 2025. "A wave-resolving two-dimensional vertical Lagrangian approach to model microplastic transport in nearshore waters based on TrackMPD 3.0", Geosci. Model Dev., 18, 319–336, https://doi.org/10.5194/gmd-18-319-2025 

The TrackMPD code is licensed under GPL (GNU General Public License). In summary, this means that the code is open source and may be used freely for non-commercial and commercial purposes. Any alterations to the TrackMPD source code or new modules must be licensed under GPL as well. See the LICENSE.md.


# Authors

TrackMPD v2.3 and v3.0: 

* Isabel Jalón-Rojas: Main developer, responsible for the model core and the behaviour, settling, washing-off, deposition-resuspension and bedload modules. 

* Vincent Marieu: Main developer, responsible for the advection module and the parallel computation. 


TrackMPD v1:

* Isabel Jalón-Rojas: Main developer, responsible for the model core, and the behaviour, settling and washing-off modules (ijalonrojas@gmail.com).
* Erick Fredj: Main developer, responsible for the advection, dispersion and output modules based on the Particle Tracking and Analysis Toolbox (PaTATO, Fredj et al., 2016) (erick.fredj@gmail.com).
* Xiao-Hua Wang: Contributor. 


# Contributing

Contributions to TrackMPD are welcome. Contact authors via email for more information.
