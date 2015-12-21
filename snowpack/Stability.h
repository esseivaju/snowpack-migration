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
#ifndef STABILITY_H
#define STABILITY_H

#include <snowpack/StabilityAlgorithms.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Constants.h>

#include <map>
#include <string>

class InstabilityData;

typedef double (*StabMemFn)(const ElementData&, const double&);
typedef bool (*StabFnShearStrength)(const double&, const double&, const mio::Date&,
                                               ElementData&, NodeData&, StabilityData&);

/**
 * @class Stability
 * @brief This class contains the stability routines for the snowpack model.
 * TODO: update description 2009-10-20 \n
 * Stability is found for each LAYER (i.e. finite element) and INTERFACE (i.e. node).
 * Subsequently, the element data contains a variable S_dr (layer stability for
 * direct action avalanches) and the nodal data contains the varialble S_i
 * (interface stability). The station data contain S_class an overall stability
 * estimation for the profile based on hand hardness, grain class and a slab
 * characterization.  At the end, the stability index approach is compared to this
 * profile characterization to check consistency/reliability. This is a first "shot"
 * and it would be a miracle if we got it right at the very beginning.
 */
class Stability {
	public:
		Stability (const SnowpackConfig& i_cfg, const bool& i_classify_profile=false);
		
		void checkStability(const CurrentMeteo& Mdata, SnowStation& Xdata);

		static const double psi_ref, max_stability, minimum_slab, ground_rough;
		static const double min_depth_ssi, skier_depth, min_thick_crust;
		static const int sh_mod, prof_classi;
		static const size_t nmax_lemon;

	private:
		void initStability(SnowStation& Xdata);
		double setStructuralStabilityIndex(const ElementData& Edata_low, const ElementData& Edata_up,
		                                   const double& Sk, InstabilityData& SIdata);
		
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static containers
		static std::map<std::string, StabMemFn> mapHandHardness;
		static std::map<std::string, StabFnShearStrength> mapShearStrength;

		std::string strength_model, hardness_parameterization;
		double hoar_density_buried;
		bool plastic;
		bool classify_profile;
};

/**
 * @class InstabilityData
 * @brief double and integer values to pinpoint structural instability
 */
class InstabilityData {
	public:
		InstabilityData() : n_lemon(0), dhard(0.0), dgsz(0.0), ssi(Stability::max_stability) {}

		size_t n_lemon;  ///< Number of "lemons" found
		double dhard;    ///< Difference in hardness
		double dgsz;     ///< Difference in grain size
		double ssi;      ///< Sk38 + structural instabilities (dhard & dgsz)
};

#endif //End of Stability.h
