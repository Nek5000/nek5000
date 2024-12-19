Reads an exodus II (.exo) file and generates a .re2 file containing the grid,
curved sides and sideset IDs. 

* Requires the 3rd-party Exodus library and NetCDF 
* Supported element types are HEX20, TET4+WEDGE6, TET4+HEX8+WEDGE6, TET10+WEDGE15+HEX20
* Periodicity should be set during exo2nek conversion. Only support translational periodicity now.
* Conjugate Heat Transfer mesh could be constructed with separate fluid and solid exo files with conformal interface
* Sideset IDs are stored in the 5th argument of the fluid bc array in re2
* The real BC's have to be specified in the .usr file using the sideset IDs
