# Simplified numerical simulation of vibrocompaction in sand 
This is an example of the Finite Element Model used in the manuscript submitted to Transportation Geotechnics by Carlos Grandas and Merita Tafili, 2024.

The example simulates the densification of a 3D sand field by vibrators ordered in a grid layout at a distance of 4 m. The material model is a combination of a hypoplastic and a HCA model. This constitutive model is implemented in a Fortran subroutine called `UMAT`, which is compatible with Abaqus.

## License
This code is licensed under the GNU General Public License, version 3 (GPL-3.0). See the LICENSE file for more details.

##  niemunis_tools_nano: A Tensor Operations Library  
This model relies on a proprietary tensor operations library developed by Andrzej Niemuni at KIT, Karlsruhe, Germany. To use this library, please contact Andrzej.Niemunis@kit.edu for access and licensing information.

## Disclaimer 
This code is provided as-is, without warranty of any kind. While every effort has been made to ensure the correctness of the subroutines, we cannot guarantee that they are free from errors. Users should thoroughly test the code in their own environments. We are not responsible for any issues arising from the use of these subroutines.


