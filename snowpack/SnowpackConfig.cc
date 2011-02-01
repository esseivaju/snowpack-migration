/*
 *  SNOWPACK VERSION 9.x
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

#include <snowpack/SnowpackConfig.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
vector<string> SnowpackConfig::variants;
map<string,string> SnowpackConfig::defaultConfig;
const bool SnowpackConfig::__init = SnowpackConfig::initStaticData();

bool SnowpackConfig::initStaticData()
{
	variants.push_back("DEFAULT");
	variants.push_back("ANTARCTICA");
	variants.push_back("JAPAN");
	variants.push_back("CALIBRATION");

	defaultConfig["ALPINE3D"] = "0";
	defaultConfig["BACKUP_DAYS_BETWEEN"] = "365";
	defaultConfig["CHECK_SW_MODE"] = "0";
	defaultConfig["CUMSUM_MASS"] = "0";
	defaultConfig["DOORSCHOT"] = "0";
	defaultConfig["ERROR_FILE"] = "./snowpack.err"; // deprecated??
	defaultConfig["FIRST_BACKUP"] = "400.";
	defaultConfig["FIXED_HEIGHTS"] = "5";
	defaultConfig["FIXED_HN_DENSITY"] = "100.";
	defaultConfig["FIXED_RATES"] = "0";
	defaultConfig["FORCE_RH_WATER"] = "1";
	defaultConfig["HARDNESS_MODEL"] = "default";
	defaultConfig["HEIGHT_NEW_ELEM"] = "0.02";
	defaultConfig["HOAR_DENSITY_BURIED"] = "125.0";
	defaultConfig["HOAR_DENSITY_SURF"] = "100.0";
	defaultConfig["HOAR_MIN_SIZE_BURIED"] = "2.0";
	defaultConfig["HOAR_MIN_SIZE_SURF"] = "0.5";
	defaultConfig["HOAR_THRESH_RH"] = "0.97";
	defaultConfig["HOAR_THRESH_VW"] = "3.5";
	defaultConfig["JAM"] = "0";
	defaultConfig["JOIN_ELEMENTS"] = "1";
	defaultConfig["MASS_BALANCE"] = "0";
	defaultConfig["MAX_NUMBER_SENSORS"] = "5";
	defaultConfig["MAX_N_SOLUTES"] = "10";
	defaultConfig["METAMORPHISM_MODEL"] = "default";
	defaultConfig["MIN_DEPTH_SUBSURF"] = "0.07";
	defaultConfig["MULTISTREAM"] = "1";
	defaultConfig["NEW_SNOW_GRAIN_RAD"] = "0.15";
	defaultConfig["NUMBER_SLOPES"] = "1";
	defaultConfig["OUT_CANOPY"] = "0";
	defaultConfig["OUT_HAZ"] = "1";
	defaultConfig["OUT_HEAT"] = "1";
	defaultConfig["OUT_LOAD"] = "0";
	defaultConfig["OUT_LW"] = "1";
	defaultConfig["OUT_MASS"] = "1";
	defaultConfig["OUT_METEO"] = "1";
	defaultConfig["OUT_STAB"] = "1";
	defaultConfig["OUT_SW"] = "1";
	defaultConfig["OUT_T"] = "1";
	defaultConfig["PERP_TO_SLOPE"] = "0";
	defaultConfig["PLASTIC"] = "0";
	defaultConfig["PRECIP_RATES"] = "1";
	defaultConfig["RESEARCH"] = "1";
	defaultConfig["STRENGTH_MODEL"] = "default";
	defaultConfig["SURFACECODE"] = "NEUMANN_BC";
	defaultConfig["THRESH_RAIN"] = "1.2";
	defaultConfig["THRESH_RH"] = "0.5";
	defaultConfig["T_CRAZY_MAX"] = "340.";
	defaultConfig["T_CRAZY_MIN"] = "210.";
	defaultConfig["VARIANT"] = variants.at(0);
	defaultConfig["VISCOSITY_MODEL"] = "default";
	defaultConfig["WET_LAYER"] = "1";
	defaultConfig["PREVAILING_WIND_DIR"] = "0.";
	defaultConfig["WIND_SCALING_FACTOR"] = "1.0";

	return true;
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

SnowpackConfig::~SnowpackConfig() {}

SnowpackConfig::SnowpackConfig(const std::string& i_filename) : Config(i_filename)
{
	string variant = get("VARIANT", "Parameters", Config::nothrow);
	int enforce_measured_snow_heights = get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Parameters");

	addKey("MINIMUM_L_ELEMENT", "Parameters", "0.0025"); //Minimum element length (m)
	double minimum_l_element = get("MINIMUM_L_ELEMENT", "Parameters"); 

	if (enforce_measured_snow_heights) {
		addKey("HEIGHT_NEW_ELEM", "Parameters", "0.02");
	} else {
		stringstream ss;
		double tmp = 2.0 * minimum_l_element;
		ss << tmp;
		addKey("HEIGHT_NEW_ELEM", "Parameters", ss.str());
	}

	if ((variant == "") || (variant == "DEFAULT")){

	} else if (variant == "JAPAN"){

		string strength_model = get("STRENGTH_MODEL", "Parameters", Config::nothrow); 
		if (strength_model == "") addKey("STRENGTH_MODEL", "Parameters", "NIED");

		string metamorphism_model = get("METAMORPHISM_MODEL", "Parameters", Config::nothrow); 
		if (metamorphism_model == "") addKey("METAMORPHISM_MODEL", "Parameters", "NIED");

		string viscosity_model = get("VISCOSITY_MODEL", "Parameters", Config::nothrow); 
		if (viscosity_model == "") addKey("VISCOSITY_MODEL", "Parameters", "VS_KOJIMA");

	} else if (variant == "ANTARCTICA"){

		addKey("MINIMUM_L_ELEMENT", "Parameters", "0.0001"); //Minimum element length (m)
		minimum_l_element = get("MINIMUM_L_ELEMENT", "Parameters");

		string hoar_density_buried = get("HOAR_DENSITY_BURIED", "Parameters", Config::nothrow);
		if (hoar_density_buried == "") addKey("HOAR_DENSITY_BURIED", "Parameters", "200.0");

		string viscosity_model = get("VISCOSITY_MODEL", "Parameters", Config::nothrow); 
		if (viscosity_model == "") addKey("VISCOSITY_MODEL", "Parameters", "VS_CALIBRATION");

		string fixed_hn_density = get("FIXED_HN_DENSITY", "Parameters", Config::nothrow);
		if (fixed_hn_density == "") addKey("FIXED_HN_DENSITY", "Parameters", "300.");

		string force_rh_water = get("FORCE_RH_WATER", "Parameters", Config::nothrow); 
		if (force_rh_water == "") addKey("FORCE_RH_WATER", "Parameters", "0");

		string thresh_rh = get("THRESH_RH", "Parameters", Config::nothrow); 
		if (thresh_rh == "") addKey("THRESH_RH", "Parameters", "0.7");

		if ( !enforce_measured_snow_heights) {
			stringstream ss;
			double tmp = 1.1 * minimum_l_element;
			ss << tmp;
			addKey("HEIGHT_NEW_ELEM", "Parameters", ss.str());
		}
		addKey("FIRST_BACKUP", "Parameters", "1200.");

		addKey("FIXED_HEIGHTS", "Parameters", "7");
		addKey("FIXED_RATES", "Parameters", "7");
		addKey("MAX_NUMBER_SENSORS", "Parameters", "23");
		addKey("MIN_DEPTH_SUBSURF", "Parameters", "0.0");
		addKey("T_CRAZY_MIN", "Parameters", "165.");
		addKey("T_CRAZY_MAX", "Parameters", "300.");
		addKey("NEW_SNOW_GRAIN_RAD", "Parameters", "0.1");

	} else if (variant == "CALIBRATION"){

		string viscosity_model = get("VISCOSITY_MODEL", "Parameters", Config::nothrow); 
		if (viscosity_model == "") addKey("VISCOSITY_MODEL", "Parameters", "VS_CALIBRATION");

		addKey("FIXED_HEIGHTS", "Parameters", "5");
		addKey("FIXED_RATES", "Parameters", "0");
		addKey("MAX_NUMBER_SENSORS", "Parameters", "8");
		addKey("MIN_DEPTH_SUBSURF", "Parameters", "0.0");

	} else {
		throw UnknownValueException("Unknown variant " + variant, AT);
	}

	//For all parameters not set by the user or by the initialization above, the default values apply
	for(map<string,string>::const_iterator it = defaultConfig.begin(); it != defaultConfig.end(); it++){
		//Loop through defaultConfig and check whether user has set the parameter, if not add default value
		string value = get(it->first, "Parameters", Config::nothrow);
		if (value == "") addKey(it->first, "Parameters", it->second);
	}

	/**
	 * @name Defines how energy and mass balance are computed \n
	 * - AVGSUM_TIME_SERIES == 1 \n
	 *   Energy and mass fluxes are averaged and cumulated over TS_DAYS_BETWEEN, respectively. \n
	 *   Otherwise, instantaneous energy fluxes and mass fluxes cumulated over the last computation
	 *   time step (CALCULATION_STEP_LENGTH) are dumped. \n
	 *   NOTE: Precipitations and Erosion are always given in rates per TS_DAYS_BETWEEN interval (h-1)
	 * - CUMSUM_MASS == 1 \n
	 *   Mass fluxes are cumulated over whole run period.
	 * - WARNING: In operational mode and if NUMBER_SLOPES > 1, the above two values are always unset!
	 */
	//@{
	int nSlopes = get("NUMBER_SLOPES", "Parameters");
	if (nSlopes < 2){
		addKey("AVGSUM_TIME_SERIES", "Parameters", "1");
	} else {
		addKey("AVGSUM_TIME_SERIES", "Parameters", "0");
	} 

	/**
	 * @brief Hazard data interval in units of CALCULATION_STEP_LENGTH \n
	 * WARNING: In operational mode, this has to result in a 30 min interval!
	 * It is a matter of consitency. If you change this, a big mess will result!!!
	 */
	double calculation_step_length = get("CALCULATION_STEP_LENGTH", "Parameters");

	string hazard_steps_between = get("HAZARD_STEPS_BETWEEN", "Parameters", Config::nothrow);
	if (hazard_steps_between == "") {
		stringstream ss;
		int tmp = (int)(30./calculation_step_length + 0.5);
		ss << tmp;
		addKey("HAZARD_STEPS_BETWEEN", "Parameters", ss.str());
	}
}

void checkUserConfiguration(mio::Config& /*cfg*/)
{
	/*
	if ( !((T_INTERNAL > -1) && (T_INTERNAL <= MAX_NUMBER_SENSORS)) ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "T_INTERNAL=%d out of range (0, %d)", T_INTERNAL, MAX_NUMBER_SENSORS);
		return ERROR;
	}
	NUMBER_SENSORS = FIXED_HEIGHTS + FIXED_RATES;
	if ( !((NUMBER_SENSORS > -1) && (NUMBER_SENSORS <= MAX_NUMBER_SENSORS)) ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "%d FIXED_HEIGHTS + %d FIXED_RATES out of range (0, %d)", FIXED_HEIGHTS, FIXED_RATES, MAX_NUMBER_SENSORS);
		return ERROR;
	}


	// Print some infos to stdout: research mode
	if ( !MEAS_TSS && CHANGE_BC ){
		prn_msg(__FILE__, __LINE__, "wrn",	-1., "(!MEAS_TSS && CHANGE_BC) == 1");
		prn_msg(__FILE__, __LINE__, "msg",	-1., "Using Neumann boundary conditions because no measured TSS is available");
		CHANGE_BC = 0;
	}
	if ( !SNP_SOIL && CANOPY ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "Canopy Model is used WITHOUT soil data!");
	}


	// Check a few settings - operational mode
	if ( SNOW_REDISTRIBUTION && ((NUMBER_SLOPES > 5) || ((NUMBER_SLOPES-1)%4 != 0)) ) {
	prn_msg(__FILE__, __LINE__, "err", -1., "NUMBER_SLOPES (%d) not compatible with SNOW_REDISTRIBUTION", NUMBER_SLOPES);
		exit(EXIT_FAILURE);
	}
	if ( CANOPY ) {
		CANOPY = 0;
		prn_msg(__FILE__, __LINE__, "wrn", -1., "CANOPY was set! You may have run into troubles! Reset to 0");
	}



	*/
}

