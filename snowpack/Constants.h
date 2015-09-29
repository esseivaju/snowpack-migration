/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>

/**
 * @file Constants.h
 * @version 10.02
 * This module defines constants needed for the 1d snowpack model
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

/**
 * @brief _VERSION is given as a compilation flag to tell us what is the version number \n
 * Please only use SN_VERSION in the code
 */
#ifndef SN_VERSION
	#define SN_VERSION _VERSION
#endif

#define MAX_STRING_LENGTH 256
#define MAX_LINE_LENGTH 6000

/// @brief Initial value for stability parameter
#define INIT_STABILITY 999.

namespace Constants {
	const double undefined = -999.; ///<This is the snowpack undefined value
	const int iundefined = -999; ///<This is the snowpack undefined value
	const size_t stundefined = static_cast<size_t>(-999);
	const double min_percent_values = 0.9;

	const double pi = 3.14159265358979323846;
	const double g  = 9.80665;
	const double stefan_boltzmann = 5.67051e-8;

	const double solcon = 1366.1; ///< Total Solar Irradiance (W m-2) (Froehlich, 2006)

	const double gas_constant = 461.9; ///< (J mol-1 K-1)
	const double gas_constant_air = 287.0; ///< for air (J kg-1 K-1)
	const double gas_constant_mol = 8.31;  ///< (J mol-1 K-1)

	/// @name Albedo (1)
	//@{
	const double min_albedo = 0.3;
	const double max_albedo = 0.99;
	const double glacier_albedo = 0.3;
	//@}
	/// @name Emissivity (1)
	//@{
	const double emissivity_snow = 0.98;
	//@}
	/// @name Density (kg m-3)
	//@{
	const double density_air = 1.1; ///< Approximation: use ideal gas law
	const double density_ice = 917.0; ///< At T = 0 degC
	const double density_water = 1000.0; ///<  At T = 0 degC
	const double density_vapor = 1.0; ///< Approximation: use ideal gas law
	/// @name Max and min densities (kg m-3)
	//@{
	const double min_rho = 5.0;
	const double max_rho = 1000.0;
	//@}
	//@}

	///@name Specific heat (J kg-1 K-1)
	//@{
	const double specific_heat_ice = 2100.0; ///< at T = 0 degC
	const double specific_heat_water = 4190.0; ///< at T = 0 degC
	const double specific_heat_air = 1004.67; ///< see Stull "Meteorology for scientists and engineers" p44
	//@}

	///@name Thermal conductivity of snow components ice, water and air (W m-1 K-1)
	//@{
	const double conductivity_ice = 2.2; ///< (W m-1 K-1)
	const double conductivity_water = 0.598; ///< (W m-1 K-1)
	const double conductivity_air = 0.026; ///< (W m-1 K-1)
	//@}

	///@name Vapor Diffusion in Air and Snow (m2 s-1)
	//@{
	///@brief This value was taken from: Colbeck, S.C., 1993. The vapor diffusion coefficient for snow, Water Resources Research, 29(1)
	const double diffusion_coefficient_in_air = 22.0e-6;
	const double diffusion_coefficient_in_snow = 85.0e-6; ///< It is larger (see Colbeck) according to work prior to 2008!
	//@}

	///@name Phase change constants
	//@{
	const double melting_tk = 273.15; ///< (K)
	const double freezing_tk = 273.15; ///< (K)
	const double triple_point_t = 273.16; ///< (K)
	const double triple_point_p = 611.73; ///< (Pa)
	const double lh_sublimation = 2.838e6; ///< (J kg-1) (solid to vapor)
	const double lh_vaporization = 2.504e6; ///< (J kg-1) (liquid to vapor)
	const double lh_fusion = 3.34e5; ///< (J kg-1) (solid to liquid)
	//@}

	///@name Numerical Constants
	//@{
	const double min_slope_angle = 3.; ///< Smallest angle a flat field may show (deg)

	/// @brief Small numbers for comparisons and to avoid divisions by zero
	const double eps = 1.e-6;
	const double eps2 = eps*eps;
	const double big = 1.0e200; ///< used for Dirichlet boundary conditions
	//@}
}

/**
 * @name MACRO definitions (to be replaced by proper functions some day)
 * TODO Move this definitions to Util.h, as functions to Util.c
 */
//@{
/// @name Elementary math functions
//@{
#define MAX(x,y)    (((x) < (y)) ? (y) : (x))
#define MIN(x,y)    (((x) < (y)) ? (x) : (y))
//@}
/// @name Conversion macros
//@{
/// @brief For lengths
#define M_TO_CM( l ) ( (l) * 100. )  // meter to centimeters
#define M_TO_MM( l ) ( (l) * 1000. ) // meter to millimeters
#define CM_TO_M( l ) ( (l) / 100. )  // centimeter to meter
#define MM_TO_M( l ) ( (l) / 1000. ) // millimeter to meter
#define MM_TO_CM( l ) ( (l) / 10. )  // millimeter to centimeter
/// @brief For time
#define D_TO_H( t ) ( (t) * 24. )    // day to hours
#define D_TO_M( t ) ( (t) * 1440. )  // day to minutes
#define D_TO_S( t ) ( (t) * 86400. ) // day to seconds
#define H_TO_D( t ) ( (t) / 24.0 )   // hour to day
#define H_TO_M( t ) ( (t) * 60. )    // hour to minutes
#define H_TO_S( t ) ( (t) * 3600. )  // hour to seconds
#define M_TO_D( t ) ( (t) / 1440. )  // minute to day
#define M_TO_H( t ) ( (t) / 60. )    // minute to hour
#define M_TO_S( t ) ( (t) * 60. )    // minute to seconds
#define S_TO_D( t ) ( (t) / 86400. ) // second to day
#define S_TO_H( t ) ( (t) / 3600. )  // second to hour
#define S_TO_M( t ) ( (t) / 60. )    // second to minute
//@}
//@}

#endif
