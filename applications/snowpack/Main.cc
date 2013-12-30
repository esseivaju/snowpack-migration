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
#include <snowpack/libsnowpack.h>
#include <meteoio/MeteoIO.h>

#include <iostream>
#include <string>
#include <sstream>

#ifdef _MSC_VER
	/*
	This software contains code under BSD license (namely, getopt for Visual C++).
	Therefore, this product includes software developed by the University of
	California, Berkeley and its contributors when compiling with Visual C++.
	*/
	#include "getopt.h"
#else
	//#include <unistd.h> //for getopt
	#include <getopt.h> //for getopt_long
#endif

using namespace std;
using namespace mio;

#ifdef DEBUG_ARITHM
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif
	#ifndef __USE_GNU
		#define __USE_GNU
	#endif
	#include <fenv.h>
#endif


/**
 * @class Slope a C. Fierz class ;-)
 */
class Slope {

	public:
		Slope(const mio::Config& cfg);

		double prevailing_wind_dir;
		unsigned int nSlopes;
		unsigned int mainStation;  ///< main station, flat field or slope
		unsigned int sector;       ///< main station (0) or current slope sector (1:nSlopes)
		unsigned int first;        ///< first virtual slope station in computing sequence
		unsigned int luv;
		unsigned int lee;
		bool north, south;
		bool snow_erosion, mainStationDriftIndex;
		bool snow_redistribution, luvDriftIndex;

		int getSectorDir(const double& dir_or_expo) const;
		void setSlope(const unsigned int slope_sequence, vector<SnowStation>& vecXdata, double& wind_dir);

	private:
		double sector_width;       ///< width of slope sector: 360./MAX(1, nSlopes-1) deg
};

/**
 * @class Cumsum a C. Fierz class ;-)
 * To cumulate various mass fluxes over either hazard or time series time steps
 */
class Cumsum {
	
	public:
		Cumsum(const unsigned int nSlopes);

		double precip;
		double drift, snow, runoff, rain;
		vector<double> erosion; // Cumulated eroded mass; dumped to file as rate
};

/************************************************************
 * static section                                           *
 ************************************************************/

//Global variables in this file:
string cfgfile = "io.ini";
string mode = "RESEARCH";
mio::Date dateEnd;
vector<string> vecStationIDs;

/// @brief Main control parameters
struct MainControl
{
	double Duration;     ///< Duration of run (s)
	size_t nStep;        ///< Time step number
	size_t nAvg;         ///< Number of calculation time steps to average fluxes etc.
	size_t HzStep;       ///< Hazard step number (should be half of nStep in operational mode)
	bool   TsDump;       ///< Flag for time series dump
	bool   HzDump;       ///< Calculation of hazard information will be performed
	bool   PrDump;       ///< Flag for profile dump
	bool   XdataDump;    ///< Backup of Xdata will be performed
	bool   sdbDump;      ///< Dump to data base if required in operational mode
	bool   resFirstDump; ///< Flag to dump initial state of snowpack
};

/************************************************************
 * non-static section                                       *
 ************************************************************/

Slope::Slope(const mio::Config& cfg)
       : prevailing_wind_dir(0.), nSlopes(0), mainStation(0), sector(0),
         first(1), luv(0), lee(0),
         north(false), south(false),
         snow_erosion(false), mainStationDriftIndex(false),
         snow_redistribution(false), luvDriftIndex(false),
         sector_width(0)
{
	cfg.getValue("NUMBER_SLOPES", "SnowpackAdvanced", nSlopes);
	cfg.getValue("SNOW_EROSION", "SnowpackAdvanced", snow_erosion);
	stringstream ss;
	ss.str("");
	ss << nSlopes;
	cfg.getValue("SNOW_REDISTRIBUTION", "SnowpackAdvanced", snow_redistribution);
		if (snow_redistribution && !(nSlopes > 1 && nSlopes % 2 == 1))
			throw mio::IOException("Please set NUMBER_SLOPES to 3, 5, 7, or 9 with SNOW_REDISTRIBUTION set! (nSlopes="+ss.str()+")", AT);
	cfg.getValue("PREVAILING_WIND_DIR", "SnowpackAdvanced", prevailing_wind_dir, mio::IOUtils::nothrow);
	sector_width = 360. / static_cast<double>(MAX(1, nSlopes-1));
}

/**
 * @brief Determine either direction of blowing wind or slope exposition\n
 *        NOTE that station slope.first always corresponds to the prevailing wind direction
 * @param dir_or_expo direction of wind or exposition
 **/
int Slope::getSectorDir(const double& dir_or_expo) const
{
	double dir = dir_or_expo;
	if (dir > 360.) dir -= 360.;
	else if (dir < 0.) dir += 360.;
	const unsigned int sectorDir = (unsigned int)((floor((dir + 0.5*sector_width)/sector_width)) + 1);
	if (sectorDir >= nSlopes) return 1;
	else return sectorDir;
}

/**
 * @brief Set slope variables
 * @param slope_sequence computation sequence for slopes
 * @param vecXdata
 * @param wind_dir direction of wind
 **/
void Slope::setSlope(const unsigned int slope_sequence, vector<SnowStation>& vecXdata, double& wind_dir)
{
	mainStationDriftIndex = false;
	luvDriftIndex = false;
	switch (slope_sequence) {
	case 0:
		for (size_t kk=0; kk<nSlopes; kk++) {
			vecXdata[kk].windward = false;
			vecXdata[kk].rho_hn   = 0.;
			vecXdata[kk].hn       = 0.;
		}
		if (nSlopes > 1) {
			luv = getSectorDir(wind_dir - prevailing_wind_dir);
			vecXdata[luv].windward = true;
			lee = (luv + nSlopes/2) % (nSlopes-1);
			if (lee == 0) lee = nSlopes - 1;
		} else {
			//requesting slope 0 of 0 expositions
			luv = lee = 0;
		}
		sector = mainStation;
		mainStationDriftIndex = ((nSlopes == 1) && snow_erosion);
		break;
	case 1:
		sector = luv;
		luvDriftIndex = snow_redistribution;
		break;
	default:
		sector++;
		if (sector == nSlopes) sector = 1;
	}
	north = (vecXdata[sector].meta.getSlopeAngle() > 0. && vecXdata[sector].meta.getAzimuth() == 0.);
	south = (vecXdata[sector].meta.getSlopeAngle() > 0. && vecXdata[sector].meta.getAzimuth() == 180.);
}

Cumsum::Cumsum(const unsigned int nSlopes)
        : precip(0.),
          drift(0.), snow(0.), runoff(0.), rain(0.),
          erosion(nSlopes, 0.)
{}

void Usage(const string& programname)
{
#ifdef _MSC_VER
	cout << "This version of Snowpack uses a BSD-licensed port of getopt for Visual C++. " << endl
		<< "It therefore includes software developed by the University of "
		<< "California, Berkeley and its contributors." << endl;
#endif
	cout << "Snowpack version " << _VERSION << " compiled on " << __DATE__ << " " << __TIME__ << endl
		<< "\tLibsnowpack " << snowpack::getLibVersion() << endl
		<< "\tMeteoIO " << mio::getLibVersion() << endl;

	cout << "Usage: " << programname << endl
		<< "\t-e, --enddate=YYYY-MM-DDTHH:MM (e.g.:2008-08-11T09:00)" << endl
		<< "\t[-c, --config=<ini file> (e.g. io.ini)]" << endl
		<< "\t[-m, --mode=[operational or research] (default: research)]" << endl
		<< "\t[-s, --stations=[comma delimited stationnames] (e.g. DAV2,WFJ2)] (NOTE: ONLY in operational mode)" << endl
		<< "\t[-h, --help] Print help message and version information]" << endl << endl;

	cout << "Example: " << programname << " -c io.ini -m research -e 1996-06-17T00:00" << endl << endl;
	exit(1);
}

void parseCmdLine(int argc, char **argv, string& end_date_str)
{
	int longindex=0, opt=-1;
	bool setEnd = false;

	struct option long_options[] =
	{
		{"enddate", required_argument, 0, 'e'},
		{"mode", required_argument, 0, 'm'},
		{"config", required_argument, 0, 'c'},
		{"stations", required_argument, 0, 's'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while ((opt=getopt_long( argc, argv, ":e:m:c:s:h", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0:
			break;
		case 'e': {
			end_date_str=string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		}
		case 'm':
			mode = string(optarg);
			mio::IOUtils::toUpper(mode);
			break;
		case 'c':
			cfgfile = string(optarg);
			break;
		case 's':
			mio::IOUtils::readLineToVec(string(optarg), vecStationIDs, ',');
			break;
		case ':': //operand missing
			cerr << endl << "[E] Command line parameter '-" << char(optopt) << "' requires an operand" << endl;
			Usage(string(argv[0]));
			break;
		case 'h':
			Usage(string(argv[0]));
			break;
		case '?':
			cerr << endl << "[E] Unknown argument detected" << endl;
			Usage(string(argv[0]));
			break;
		default:
			cerr << endl << "[E] getopt returned character code " <<  opt << endl;
			Usage(string(argv[0]));
		}
	}

	if (!setEnd) {
		cerr << endl << "[E] You must specify an enddate for the simulation!" << endl;
		Usage(string(argv[0]));
	}
}

void editMeteoData(mio::MeteoData& md, const string& variant)
{
	// Since we cannot deal with precipitation nodata, we set it to zero (HACK)
	if (md(MeteoData::HNW) == mio::IOUtils::nodata)
		md(MeteoData::HNW) = 0.0;

	if (md(MeteoData::VW) == mio::IOUtils::nodata)
		md(MeteoData::VW) = 1.0; // if no wind measurement exists assume 1 m/s; ori: 3 m/s

	if (md(MeteoData::DW) == mio::IOUtils::nodata)
		md(MeteoData::DW) = 0.;

	if (md(MeteoData::TSG) == mio::IOUtils::nodata)
		md(MeteoData::TSG) = 273.15;

	//Add the atmospheric emissivity as a parameter
	if (!md.param_exists("EA")) md.addParameter("EA");
		md("EA") = SnLaws::AirEmissivity(md, variant);

	// Snow stations without separate wind station use their own wind for local drifting and blowing snow
	if (!md.param_exists("VW_DRIFT")) {
		md.addParameter("VW_DRIFT");
		md("VW_DRIFT") = md(MeteoData::VW);
	}
	if (!md.param_exists("DW_DRIFT")) {
		md.addParameter("DW_DRIFT");
		md("DW_DRIFT") = md(MeteoData::DW);
	}
}

// Return true if snowpack can compute the next timestep, else false
bool validMeteoData(const mio::MeteoData& md, const string& StationName, const string& variant)
{
	bool miss_ta=false, miss_rh=false, miss_precip=false, miss_rad=false;
	bool miss_ea=false;

	if (md(MeteoData::TA) == mio::IOUtils::nodata)
		miss_ta=true;
	if (md(MeteoData::RH) == mio::IOUtils::nodata)
		miss_rh=true;
	if ((variant != "ANTARCTICA")
	        && ((md(MeteoData::ISWR) == mio::IOUtils::nodata) && (md(MeteoData::RSWR) == mio::IOUtils::nodata)))
		miss_rad=true;
	if ((md(MeteoData::HNW) == mio::IOUtils::nodata) && (md(MeteoData::HS) == mio::IOUtils::nodata))
		miss_precip=true;
	if (md("EA") == mio::IOUtils::nodata)
		miss_ea=true;

	if(miss_ta || miss_rh || miss_rad || miss_precip || miss_ea) {
		mio::Date now;
		now.setFromSys();
		cerr << "[E] [" << now.toString(mio::Date::ISO) << "] ";
		cerr << StationName << " missing { ";
		if(miss_ta) cerr << "TA ";
		if(miss_rh) cerr << "RH ";
		if(miss_rad) cerr << "radiation ";
		if(miss_precip) cerr << "precipitations ";
		if(miss_ea) cerr << "ea ";
		cerr << "} on " << md.date.toString(mio::Date::ISO) << "\n";
		return false;
	}
	return true;
}

void copyMeteoData(const mio::MeteoData& md, CurrentMeteo& Mdata,
                   const double prevailing_wind_dir, const double wind_scaling_factor)
{
	Mdata.date   = Date::rnd(md.date, 1);
	Mdata.ta     = md(MeteoData::TA);
	Mdata.rh     = md(MeteoData::RH);
	if (md.param_exists("RH_AVG"))
		Mdata.rh_avg = md("RH_AVG");
	Mdata.vw     = md(MeteoData::VW);
	Mdata.dw     = md(MeteoData::DW);
	Mdata.vw_max = md(MeteoData::VW_MAX);
	if (md.param_exists("VW_AVG"))
		Mdata.vw_avg = md("VW_AVG");

	Mdata.vw_drift = md("VW_DRIFT");
	if (Mdata.vw_drift != mio::IOUtils::nodata) Mdata.vw_drift *= wind_scaling_factor;
	Mdata.dw_drift = md("DW_DRIFT");
	if (Mdata.dw_drift == mio::IOUtils::nodata) Mdata.dw_drift = prevailing_wind_dir;

	Mdata.iswr   = md(MeteoData::ISWR);
	Mdata.rswr   = md(MeteoData::RSWR);

	Mdata.ea  = md("EA");
	Mdata.tss = md(MeteoData::TSS);
	if (md.param_exists("TSS_A12H") && (md("TSS_A12H") != mio::IOUtils::nodata))
		Mdata.tss_a12h = md("TSS_A12H");
	else
		Mdata.tss_a12h = Constants::undefined;
	if (md.param_exists("TSS_A24H") && (md("TSS_A24H") != mio::IOUtils::nodata))
		Mdata.tss_a24h = md("TSS_A24H");
	else
		Mdata.tss_a24h = Constants::undefined;
	Mdata.ts0 = md(MeteoData::TSG);
	Mdata.hnw = md(MeteoData::HNW);
	Mdata.hs = md(MeteoData::HS);
	if (md.param_exists("HS_A3H") && (md("HS_A3H") != mio::IOUtils::nodata))
		Mdata.hs_a3h = md("HS_A3H");
	else
		Mdata.hs_a3h = Constants::undefined;

	// Add measured new snow density if available
	if (md.param_exists("RHO_HN"))
		Mdata.rho_hn = md("RHO_HN");

	// Add advective heat (for permafrost) if available
	if(md.param_exists("ADV_HEAT"))
		Mdata.adv_heat = md("ADV_HEAT");
}

/**
 * @brief Make sure that both short wave fluxes get at least a "realistic" value but measured albedo only if both fluxes are measured
 * @note To be done only for flat field or single slope station
 * @param Mdata
 * @param Xdata
 * @param slope
 */
void setShortWave(CurrentMeteo& Mdata, const SnowStation& Xdata, const bool& iswr_is_net)
{
	if ((Mdata.iswr > 5.) && (Mdata.rswr > 3.) && !iswr_is_net)
		Mdata.mAlbedo = Mdata.rswr / Mdata.iswr;
	else
		Mdata.mAlbedo = Constants::undefined;

	const double cAlbedo = Xdata.Albedo;
	
	if(iswr_is_net) {
		const double netSW = Mdata.iswr;
		if(netSW==0.) { //this should only happen at night
			Mdata.iswr = 0.;
			Mdata.rswr = 0.;
			return;
		}
		Mdata.iswr = netSW / (1. - cAlbedo);
		Mdata.rswr = netSW / (1./cAlbedo - 1.);
		return;
	}
	
	if (Mdata.iswr == mio::IOUtils::nodata)
		Mdata.iswr = Mdata.rswr / Xdata.Albedo;
	if (Mdata.rswr == mio::IOUtils::nodata)
		Mdata.rswr = Mdata.iswr * Xdata.Albedo;
}

//for a given config (that can be altered) and original meteo data, prepare the snowpack data structures
//This means that all tweaking of config MUST be reflected in the config object
void dataForCurrentTimeStep(CurrentMeteo& Mdata, SurfaceFluxes& surfFluxes, vector<SnowStation>& vecXdata,
                            const Slope& slope, SnowpackConfig& cfg,
                            SunObject &sun,
                            double& precip, const double& lw_in, const double hs_a3hl6,
                            double& tot_mass_in)
{
	SnowStation &currentSector = vecXdata[slope.sector]; //alias: the current station
	const bool isMainStation = (slope.sector == slope.mainStation);
	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	const bool perp_to_slope = cfg.get("PERP_TO_SLOPE", "SnowpackAdvanced");
	const bool iswr_is_net = cfg.get("ISWR_IS_NET", "Input");
	if (Mdata.tss == mio::IOUtils::nodata) {
		// NOTE In case CHANGE_BC is set, this leads to degraded computation, that is, use parameterized
		//      incoming long wave with NEUMANN BC; it's better than nothing if no TSS is available!
		cfg.addKey("MEAS_TSS", "Snowpack", "false");
	}

	Meteo meteo(cfg);

	// Reset Surface and Canopy Data to zero if you seek current values
	const bool avgsum_time_series = cfg.get("AVGSUM_TIME_SERIES", "Output");
	if (!avgsum_time_series) {
		const bool cumsum_mass = cfg.get("CUMSUM_MASS", "Output");
		surfFluxes.reset(cumsum_mass);
		if (useCanopyModel)
			currentSector.Cdata.reset(cumsum_mass);

		const bool mass_balance = cfg.get("MASS_BALANCE", "SnowpackAdvanced");
		if (mass_balance) {
			// Do an initial mass balance check
			if (!massBalanceCheck(currentSector, surfFluxes, tot_mass_in))
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Mass error during initial check!");
		}
	}
	// Reset surfFluxes.drift and surfFluxes.mass[MS_WIND] anyway since you use these to transport eroded snow
	surfFluxes.drift = 0.;
	surfFluxes.mass[SurfaceFluxes::MS_WIND] = 0.;

	if (isMainStation) {
		// Check for growing grass
		if (!meteo.compHSrate(Mdata, currentSector, hs_a3hl6))
			cfg.addKey("DETECT_GRASS", "SnowpackAdvanced", "false");

		// Set iswr/rswr and measured albedo
		setShortWave(Mdata, currentSector, iswr_is_net);
		meteo.compRadiation(currentSector, sun, cfg, Mdata);
	} else { // Virtual slope
		cfg.addKey("CHANGE_BC", "Snowpack", "false");
		cfg.addKey("MEAS_TSS", "Snowpack", "false");
		Mdata.tss = Constants::undefined;
		cfg.addKey("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "true");
		cfg.addKey("DETECT_GRASS", "SnowpackAdvanced", "false");
	}

	const int sw_mode = static_cast<int>(cfg.get("SW_MODE", "Snowpack")) % 10; //it must be after calling compRadiation!

	// Project irradiance on slope; take care of measured snow depth and/or precipitations too
	if (!perp_to_slope) {
		meteo.radiationOnSlope(currentSector, sun, Mdata, surfFluxes);
		if ( ((sw_mode == 1) || (sw_mode == 2))
			&& (currentSector.meta.getSlopeAngle() > Constants::min_slope_angle)) { // Do not trust blindly measured RSWR on slopes
			cfg.addKey("SW_MODE", "Snowpack", "0"); // as Mdata.iswr is the sum of dir_slope and diff
		}
		meteo.projectPrecipitations(currentSector.meta.getSlopeAngle(), Mdata.hnw, Mdata.hs);
	}

	// Find the Wind Profile Parameters, w/ or w/o canopy; take care of canopy
	meteo.compMeteo(Mdata, currentSector);

	if (isMainStation) {
		// Update precipitation memory of main station
		precip += Mdata.hnw;
		if (Mdata.hs != mio::IOUtils::nodata) {
			currentSector.mH = Mdata.hs + currentSector.Ground;
		}
	} else { // Virtual slope
		currentSector.mH = Constants::undefined;
		
		// A) Compute depth of snowfall (hn*) and new snow density (rho_hn*)
		double hn_slope = 0., rho_hn_slope = SnLaws::min_hn_density;
		if (vecXdata[slope.mainStation].hn > 0.) {
			// Assign new snow depth and density from station field (usually flat)
			hn_slope = vecXdata[slope.mainStation].hn * currentSector.cos_sl;
			rho_hn_slope = vecXdata[slope.mainStation].rho_hn;
		}
		/*
		 * Snow redistribution on slopes: Add windward eroded snow to lee slope
		 * These are very important lines: Note that deposition is treated here (lee)
		 * while erosion is treated in SnowDrift.c (windward).
		*/
		if (slope.snow_redistribution && (slope.sector == slope.lee)) {
			// If it is not snowing, use surface snow density on windward slope
			if (!(hn_slope > 0.)) {
				rho_hn_slope = vecXdata[slope.luv].rho_hn;
			}
			// Add eroded mass from windward slope
			if (rho_hn_slope != 0.) {
				hn_slope += vecXdata[slope.luv].ErosionMass / rho_hn_slope;
			}
			vecXdata[slope.luv].ErosionMass = 0.;
		}
		// Update depth of snowfall on slopes.
		// This may include contributions from drifting snow eroded on the windward (luv) slope.
		if ((hn_slope > 0.) && (vecXdata[slope.mainStation].cH > 0.01)) {
			currentSector.hn = hn_slope;
			currentSector.rho_hn = rho_hn_slope;
		}

		// B) Check whether to use incoming longwave as estimated from station field
		const bool meas_incoming_longwave = cfg.get("MEAS_INCOMING_LONGWAVE", "SnowpackAdvanced");
		if (!meas_incoming_longwave) {
			double ea = 0.;
			if ((lw_in > 50.) && (lw_in < 300.)) {
				ea = Atmosphere::blkBody_Emissivity(lw_in, Mdata.ta);
			}
			if ((ea > 0.55) && (ea < 1.0)) {
				Mdata.ea = ea;
			}
		}
	}
}

/**
 * @brief determine which outputs need to be done for the current time step
 * @param mn_ctrl timestep control structure
 * @param step current time integration step
 * @param sno_step current step in the sno files (current sno profile)
 */

void getOutputControl(MainControl& mn_ctrl, const mio::Date& step, const mio::Date& sno_step,
                      const double& calculation_step_length,
                      const double& tsstart, const double& tsdaysbetween,
                      const double& profstart, const double& profdaysbetween,
                      const double& first_backup, const double& backup_days_between)
{
//HACK: put all tsstart, tsdaysbetween, etc in MainControl as well as current timestep
	const double Dstep = step.getJulian();
	const double Dsno_step = sno_step.getJulian();
	if (mn_ctrl.resFirstDump) {
		mn_ctrl.HzDump = false;
		mn_ctrl.TsDump = true;
		mn_ctrl.PrDump = true;
		mn_ctrl.resFirstDump = false;
	} else {
		// Hazard data, every half-hour
		mn_ctrl.HzDump = booleanTime(Dstep, 0.5/24., 0.0, calculation_step_length);
		// Time series (*.met)
		double bool_start = H_TO_D(tsstart);
		if ( bool_start > 0. )
			bool_start += Dsno_step;
		mn_ctrl.TsDump = booleanTime(Dstep, tsdaysbetween, bool_start, calculation_step_length);
		// Profile (*.pro)
		bool_start = H_TO_D(profstart);
		if ( bool_start > 0. )
			bool_start += Dsno_step;
		mn_ctrl.PrDump = booleanTime(Dstep, profdaysbetween, bool_start, calculation_step_length);

	}

	// Additional Xdata backup (*.<JulianDate>sno)
	const double bool_start = Dsno_step + first_backup;
	mn_ctrl.XdataDump = booleanTime(Dstep, backup_days_between, bool_start, calculation_step_length);
}

bool readSlopeMeta(mio::IOManager& io, SnowpackIO& snowpackio, SnowpackConfig& cfg, const size_t& i_stn,
                   Slope& slope, mio::Date &current_date, vector<SN_SNOWSOIL_DATA> &vecSSdata,
                   vector<SnowStation> &vecXdata, ZwischenData &sn_Zdata, CurrentMeteo& Mdata,
                   double &wind_scaling_factor, double &time_count_deltaHS)
{
	string snowfile("");
	stringstream ss;
	ss.str("");
	ss << "SNOWFILE" << i_stn+1;
	cfg.getValue(ss.str(), "Input", snowfile, mio::IOUtils::nothrow);

	//Read SSdata for every "slope" referred to as sector where sector 0 corresponds to the main station
	for (size_t sector=slope.mainStation; sector<slope.nSlopes; sector++) {
		try {
			if (sector == slope.mainStation) {
				if (snowfile == "") {
					snowfile = vecStationIDs[i_stn];
				} else {
					const size_t pos_dot = snowfile.rfind(".");
					const size_t pos_slash = snowfile.rfind("/");
					if (((pos_dot != string::npos) && (pos_dot > pos_slash)) ||
						((pos_dot != string::npos) && (pos_slash == string::npos))) //so that the dot is not in a directory name
						snowfile.erase(pos_dot, snowfile.size()-pos_dot);
				}
				snowpackio.readSnowCover(snowfile, vecStationIDs[i_stn], vecSSdata[slope.mainStation], sn_Zdata);
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Reading snow cover data for station %s",
				        vecStationIDs[i_stn].c_str());
				// NOTE (Is it a HACK?) Reading station meta data provided in meteo data and prebuffering those data
				vector<mio::MeteoData> vectmpmd;
				current_date = Date::rnd(vecSSdata[slope.mainStation].profileDate, 1);
				io.getMeteoData(current_date, vectmpmd);
				if (vectmpmd.empty())
					throw mio::IOException("No data found for station " + vecStationIDs[i_stn] + " on "
					                       + current_date.toString(mio::Date::ISO), AT);
				Mdata.setMeasTempParameters(vectmpmd[i_stn]);
				vecSSdata[slope.mainStation].meta = mio::StationData::merge(vectmpmd[i_stn].meta,
				                                vecSSdata[slope.mainStation].meta);
			} else {
				stringstream sec_snowfile;
				sec_snowfile << snowfile << sector;
				ss.str("");
				ss << vecSSdata[slope.mainStation].meta.getStationID() << sector;
				snowpackio.readSnowCover(sec_snowfile.str(), ss.str(), vecSSdata[sector], sn_Zdata);
				vecSSdata[sector].meta.position = vecSSdata[slope.mainStation].meta.getPosition();
				vecSSdata[sector].meta.stationName = vecSSdata[slope.mainStation].meta.getStationName();
			}
			vecXdata[sector].initialize(vecSSdata[sector], sector); // Generate the corresponding Xdata
		} catch (const exception& e) {
			if (sector == slope.first) {
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
				        "No virtual slopes! Computation for main station %s only!", vecStationIDs[i_stn].c_str());
				slope.nSlopes = 1;
				if ((mode == "OPERATIONAL")
					&& (vecSSdata[slope.mainStation].meta.getSlopeAngle() > Constants::min_slope_angle)) {
					cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "true");
				}
				break;
			} else {
				cout << e.what();
				throw;
			}
		}

		// Operational mode ONLY: Pass ... wind_factor, snow depth discrepancy time counter
		if ((sector == slope.mainStation) && (mode == "OPERATIONAL")) {
			wind_scaling_factor = vecSSdata[slope.mainStation].WindScalingFactor;
			time_count_deltaHS = vecSSdata[slope.mainStation].TimeCountDeltaHS;
		}
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Finished Initializing station %s", vecStationIDs[i_stn].c_str());

	//CHECK date inconsistencies between sno files
	bool dates_consistent(true);
	for (size_t sector=slope.first; sector<slope.nSlopes; sector++) {
		if (vecSSdata[sector].profileDate != vecSSdata[slope.mainStation].profileDate) {
			prn_msg(__FILE__, __LINE__, "err", mio::Date(),
				"Date virtual slope %d inconsistent with flat field %s", sector, vecStationIDs[i_stn].c_str());
			dates_consistent = false;

		}
	}
	if (!dates_consistent) return false; //go to next station

	// Do not go ahead if starting time is larger than maxtime!
	if (vecSSdata[slope.mainStation].profileDate > dateEnd) {
		prn_msg(__FILE__, __LINE__, "err", mio::Date(),
			"Starting time (%.5lf) larger than end time(%.5lf), station %s!",
			vecSSdata[slope.mainStation].profileDate.getJulian(), dateEnd.getJulian(),
			vecStationIDs[i_stn].c_str());
		return false; //goto next station
	}

	return true;
}

// SNOWPACK MAIN **************************************************************
void real_main (int argc, char *argv[])
{
	if (argc==1) Usage(string(argv[0]));

#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	const bool prn_check = false;
	mio::Timer meteoRead_timer;
	mio::Timer run_timer;
	run_timer.start();
	time_t nowSRT = time(NULL);
	MainControl mn_ctrl; //Time step control parameters

	string end_date_str("");
	parseCmdLine(argc, argv, end_date_str); //parse the command line arguments

	SnowpackConfig cfg(cfgfile);

	const double i_time_zone = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	if (end_date_str == "NOW") { //interpret user provided end date
		dateEnd.setFromSys();
		dateEnd.setTimeZone(i_time_zone);
		dateEnd.rnd(1800, mio::Date::DOWN);
	} else {
		mio::IOUtils::convertString(dateEnd, end_date_str, i_time_zone);
	}

	string outpath(""), experiment(""), db_name(""), variant("");
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant, mio::IOUtils::nothrow);

	// Add keys to perform running mean in Antarctic variant
	if (variant == "ANTARCTICA") {
		cfg.addKey("COPY::VW_AVG", "Input", "VW");
		cfg.addKey("COPY::RH_AVG", "Input", "RH");

		cfg.addKey("VW_AVG::filter1", "Filters", "mean_avg");
		cfg.addKey("VW_AVG::arg1", "Filters", "soft 101 360000");
		cfg.addKey("RH_AVG::filter1", "Filters", "mean_avg");
		cfg.addKey("RH_AVG::arg1", "Filters", "soft 101 360000");
	}

	const int tst_sw_mode = cfg.get("SW_MODE", "Snowpack"); // Test settings for SW_MODE
	if ((tst_sw_mode % 10) == 2) { //HACK: this is only for INP!
		// Make sure there is not only one of ISWR and RSWR available
		bool iswr_inp=true, rswr_inp = true;
		cfg.getValue("ISWR_INP","Input",iswr_inp,IOUtils::nothrow);
		cfg.getValue("RSWR_INP","Input",rswr_inp,IOUtils::nothrow);
		if (!(iswr_inp && rswr_inp)) {
			cerr << "[E] SW_MODE = 2: Please set both ISWR_INP and RSWR_INP to true in [Input]-section of io.ini!\n";
			exit(1);
		}
	}

	const bool useSoilLayers = cfg.get("SNP_SOIL", "Snowpack");
	bool soil_flux = false;
	if (useSoilLayers) {
		cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux);
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK w/ soil layers in %s mode", mode.c_str());
	} else {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK in %s mode", mode.c_str());
	}
	if (variant != "DEFAULT") {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Variant is '%s'", variant.c_str());
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
	        "%s compiled on %s at %s", string(argv[0]).c_str(), __DATE__, __TIME__);

	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	bool detect_grass = cfg.get("DETECT_GRASS", "SnowpackAdvanced");
	if (mode == "OPERATIONAL") {
		cfg.addKey("RESEARCH", "SnowpackAdvanced", "false");
		cfg.addKey("AVGSUM_TIME_SERIES", "Output", "false");
		cfg.getValue("DBNAME", "Output", db_name, mio::IOUtils::nothrow);
		if (useCanopyModel) {
			throw mio::IOException("Please don't set CANOPY to 1 in OPERATIONAL mode", AT);
		}
		if (!detect_grass){
			cfg.addKey("DETECT_GRASS", "SnowpackAdvanced", "true");
			detect_grass = true;
		}
	} else {
		cfg.getValue("EXPERIMENT", "Output", experiment, mio::IOUtils::nothrow);
		cfg.getValue("METEOPATH", "Output", outpath);
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Experiment : %s", experiment.c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Output dir : %s", outpath.c_str());
	}
	if (detect_grass) {
		// we need various average values of tss and hs, all for "past" windows (left)
		// Require at least one value per 3 hours
		cfg.addKey("COPY::TSS_A24H", "Input", "TSS");
		cfg.addKey("TSS_A24H::filter1", "Filters", "mean_avg");
		cfg.addKey("TSS_A24H::arg1", "Filters", "left 48 86340"); //TODO change # data required to 4

		cfg.addKey("COPY::TSS_A12H", "Input", "TSS");
		cfg.addKey("TSS_A12H::filter1", "Filters", "mean_avg");
		cfg.addKey("TSS_A12H::arg1", "Filters", "left 24 43140"); //TODO change # data required to 2

		cfg.addKey("COPY::HS_A3H", "Input", "HS");
		cfg.addKey("HS_A3H::filter1", "Filters", "mean_avg");
		cfg.addKey("HS_A3H::arg1", "Filters", "left 6 10740"); //TODO change # data required to 1
	}

	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	const double sn_dt = M_TO_S(calculation_step_length); //Calculation time step in seconds

	int nSolutes = Constants::iundefined;
	cfg.getValue("NUMBER_OF_SOLUTES", "Input", nSolutes, mio::IOUtils::nothrow);
	if (nSolutes > 0) SnowStation::number_of_solutes = static_cast<short unsigned int>(nSolutes);

	//Interval between profile backups (*.sno\<JulianDate\>) (d)
	double backup_days_between = 400.;
	cfg.getValue("BACKUP_DAYS_BETWEEN", "Output", backup_days_between, mio::IOUtils::nothrow);
	//First additional profile backup (*.sno\<JulianDate\>) since start of simulation (d)
	double first_backup = 0.;
	cfg.getValue("FIRST_BACKUP", "Output", first_backup, mio::IOUtils::nothrow);

	const bool classify_profile = cfg.get("CLASSIFY_PROFILE", "Output", mio::IOUtils::nothrow);
	const bool profwrite = cfg.get("PROF_WRITE", "Output");
	const double profstart = cfg.get("PROF_START", "Output");
	const double profdaysbetween = cfg.get("PROF_DAYS_BETWEEN", "Output");
	const bool tswrite = cfg.get("TS_WRITE", "Output");
	const double tsstart = cfg.get("TS_START", "Output");
	const double tsdaysbetween = cfg.get("TS_DAYS_BETWEEN", "Output");

	const bool precip_rates = cfg.get("PRECIP_RATES", "Output", mio::IOUtils::nothrow);
	const bool avgsum_time_series = cfg.get("AVGSUM_TIME_SERIES", "Output", mio::IOUtils::nothrow);
	const bool cumsum_mass = cfg.get("CUMSUM_MASS", "Output", mio::IOUtils::nothrow);

	//If the user provides the stationIDs - operational use case
	if (!vecStationIDs.empty()) { //This means that the user provides the station IDs on the command line
		for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
			stringstream ss;
			ss << "STATION" << i_stn+1; //For the IMIS plugin of MeteoIO, this key denotes the station id
			cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]);
		}
	}

	SnowpackIO snowpackio(cfg);

	mio::IOManager io(cfg);
	io.setMinBufferRequirements(IOUtils::nodata, 1.1); //we require the buffer to contain at least 1.1 day before the current point

	vector<StationData> accessible_stations;
	io.getStationData(dateEnd, accessible_stations); //we are retrieving meta information from MeteoIO
	if (vecStationIDs.empty()) {
		for (size_t ii=0; ii<accessible_stations.size(); ii++) {
			vecStationIDs.push_back(accessible_stations[ii].getStationID()); //HACK: accessible_stations should be directly used
		}
	}

	// START LOOP OVER ALL STATIONS
	for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
		cout << endl;
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Run on station %s", vecStationIDs[i_stn].c_str());
		run_timer.reset();
		meteoRead_timer.reset();

		Slope slope(cfg);
		Cumsum cumsum(slope.nSlopes);

		double lw_in = Constants::undefined;    // Storage for LWin from flat field energy balance

		// Used to scale wind for blowing and drifting snowpack (from statistical analysis)
		double wind_scaling_factor = cfg.get("WIND_SCALING_FACTOR", "SnowpackAdvanced");

		// Control of time window: used for adapting diverging snow depth in operational mode
		double time_count_deltaHS = 0.;

		// Snowpack data (input/output)
		vector<SN_SNOWSOIL_DATA> vecSSdata(slope.nSlopes, SN_SNOWSOIL_DATA(/*number_of_solutes*/));
		vector<SnowStation> vecXdata(slope.nSlopes, SnowStation(useCanopyModel, useSoilLayers/*, number_of_solutes*/));
		ZwischenData sn_Zdata;   // "Memory"-data, required for every operational station

		// Create meteo data object to hold interpolated current time steps
		CurrentMeteo Mdata(cfg);
		// To collect surface exchange data for output
		SurfaceFluxes surfFluxes/*(number_of_solutes)*/;
		// Boundary condition (fluxes)
		BoundCond sn_Bdata;

		mio::Date current_date;
		meteoRead_timer.start();
		if (mode == "OPERATIONAL")
			cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "false");
		const bool read_slope_status = readSlopeMeta(io, snowpackio, cfg, i_stn, slope, current_date, vecSSdata, vecXdata, sn_Zdata, Mdata, wind_scaling_factor, time_count_deltaHS);
		meteoRead_timer.stop();
		if(!read_slope_status) continue; //something went wrong, move to the next station

		memset(&mn_ctrl, 0, sizeof(MainControl));
		if (mode == "RESEARCH") {
			mn_ctrl.resFirstDump = true; //HACK to dump the initial state in research mode
			deleteOldOutputFiles(outpath, experiment, vecStationIDs[i_stn], slope.nSlopes);
			cfg.write(outpath + "/" + vecStationIDs[i_stn] + "_" + experiment + ".ini"); //output config
			current_date -= calculation_step_length/1440;
		} else {
			if (!db_name.empty() && (db_name == "sdbo" || db_name == "sdbt"))
				mn_ctrl.sdbDump = true;
		}

		SunObject sun(vecSSdata[slope.mainStation].meta.position.getLat(), vecSSdata[slope.mainStation].meta.position.getLon(), vecSSdata[slope.mainStation].meta.position.getAltitude());
		sun.setElevationThresh(0.6);
		mn_ctrl.Duration = (dateEnd.getJulian() - vecSSdata[slope.mainStation].profileDate.getJulian() + 0.5/24)*24*3600;
		vector<ProcessDat> qr_Hdata;     //Hazard data for t=0...tn
		vector<ProcessInd> qr_Hdata_ind; //Hazard data Index for t=0...tn
		Hazard hazard(cfg, mn_ctrl.Duration);
		hazard.initializeHazard(sn_Zdata.drift24, vecXdata.at(0).meta.getSlopeAngle(), qr_Hdata, qr_Hdata_ind);

		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Start simulation for %s on %s",
			vecStationIDs[i_stn].c_str(), vecSSdata[slope.mainStation].profileDate.toString(mio::Date::ISO_TZ).c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "End date specified by user: %s",
		        dateEnd.toString(mio::Date::ISO_TZ).c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Integration step length: %f min",
		        calculation_step_length);

		bool computed_one_timestep = false;
		double meteo_step_length = -1.;

		// START TIME INTEGRATION LOOP
		do {
			current_date += calculation_step_length/1440;
			mn_ctrl.nStep++;
			mn_ctrl.nAvg++;

			// Get meteo data
			vector<mio::MeteoData> vecMyMeteo;
			meteoRead_timer.start();
			io.getMeteoData(current_date, vecMyMeteo);
			if(meteo_step_length<0.) {
				std::stringstream ss2;
				meteo_step_length = io.getAvgSamplingRate();
				ss2 << meteo_step_length;
				cfg.addKey("METEO_STEP_LENGTH", "Snowpack", ss2.str());
			}
			meteoRead_timer.stop();
			editMeteoData(vecMyMeteo[i_stn], variant);
			if (!validMeteoData(vecMyMeteo[i_stn], vecStationIDs[i_stn], variant)) {
				prn_msg(__FILE__, __LINE__, "msg-", current_date, "No valid data for station %s on [%s]",
				        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str());
				current_date -= calculation_step_length/1440;
				break;
			}

			//determine which outputs will have to be done
			getOutputControl(mn_ctrl, current_date, vecSSdata[slope.mainStation].profileDate, calculation_step_length,
			                 tsstart, tsdaysbetween, profstart, profdaysbetween,
			                 first_backup, backup_days_between);
			//Radiation data
			sun.setDate(current_date.getJulian(), current_date.getTimeZone());

			std::vector<mio::MeteoData> MyMeteol3h;
			try {
				io.getMeteoData(current_date - 3.0/24.0, MyMeteol3h);  // meteo data with 3 h (left) lag
			} catch (...) {
				cerr << "[E] failed to read meteo data with 3 hours (left) lag\n";
				throw;
			}
			double hs_a3hl6;
			if (MyMeteol3h[0].param_exists("HS_A3H") && (MyMeteol3h[0]("HS_A3H") != mio::IOUtils::nodata))
				hs_a3hl6 = MyMeteol3h[0]("HS_A3H");
			else
				hs_a3hl6 = Constants::undefined;

			// START LOOP OVER ASPECTS
			for (size_t slope_sequence=0; slope_sequence<slope.nSlopes; slope_sequence++) {
				double tot_mass_in = 0.; // To check mass balance over one CALCULATION_STEP_LENGTH if MASS_BALANCE is set
				SnowpackConfig tmpcfg(cfg);
				copyMeteoData(vecMyMeteo[i_stn], Mdata, slope.prevailing_wind_dir, wind_scaling_factor);
				Mdata.copySnowTemperatures(vecMyMeteo[i_stn], slope_sequence);
				Mdata.copySolutes(vecMyMeteo[i_stn], SnowStation::number_of_solutes);
				slope.setSlope(slope_sequence, vecXdata, Mdata.dw_drift);
				dataForCurrentTimeStep(Mdata, surfFluxes, vecXdata, slope, tmpcfg,
                                       sun, cumsum.precip, lw_in, hs_a3hl6,
                                       tot_mass_in);

				// Notify user every fifteen days of date being processed
				const double notify_start = floor(vecSSdata[slope.mainStation].profileDate.getJulian()) + 15.5;
				if ((mode == "RESEARCH") && (slope.sector == slope.mainStation)
				        && booleanTime(current_date.getJulian(), 15., notify_start, calculation_step_length)) {
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					            "Station %s (%d slope(s)): advanced to %s (%f) station time",
					                vecSSdata[slope.mainStation].meta.stationID.c_str(), slope.nSlopes,
					                    current_date.toString(mio::Date::DIN).c_str(), current_date.getJulian());
				}

				// SNOWPACK model (Temperature and Settlement computations)
				Snowpack snowpack(tmpcfg); //the snowpack model to use
				Stability stability(tmpcfg, classify_profile);
				snowpack.runSnowpackModel(Mdata, vecXdata[slope.sector], cumsum.precip, sn_Bdata, surfFluxes);
				stability.checkStability(Mdata, vecXdata[slope.sector]);

				/***** OUTPUT SECTION *****/
				surfFluxes.collectSurfaceFluxes(sn_Bdata, vecXdata[slope.sector], Mdata);
				if (slope.sector == slope.mainStation) { // main station only (usually flat field)
					// Calculate consistent lw_in for virtual slopes
					if ( vecXdata[slope.mainStation].getNumberOfElements() > 0 ) {
						double k_eff, gradT;
						k_eff =
						    vecXdata[slope.mainStation].Edata[vecXdata[slope.mainStation].getNumberOfElements()-1].k[TEMPERATURE];
						gradT =
						    vecXdata[slope.mainStation].Edata[vecXdata[slope.mainStation].getNumberOfElements()-1].gradT;
						lw_in = k_eff*gradT + sn_Bdata.lw_out - sn_Bdata.qs - sn_Bdata.ql - sn_Bdata.qr;
					} else {
						lw_in = Constants::undefined;
					}
					// Deal with new snow densities
					if (vecXdata[slope.mainStation].hn > 0.) {
						surfFluxes.cRho_hn = vecXdata[slope.mainStation].rho_hn;
						surfFluxes.mRho_hn = Mdata.rho_hn;
					}
					if (slope.snow_erosion) {
						// Update drifting snow index (VI24),
						//   from erosion at the main station only if no virtual slopes are available
						if (slope.mainStationDriftIndex)
							cumulate(cumsum.drift, surfFluxes.drift);
						// Update erosion mass from main station
						// NOTE cumsum.erosion[] will be positive in case of real erosion at any time during the output time step
						if (vecXdata[slope.mainStation].ErosionMass > Constants::eps) {
							// Real erosion
							if (cumsum.erosion[slope.mainStation] > Constants::eps)
								cumsum.erosion[slope.mainStation] += vecXdata[slope.mainStation].ErosionMass;
							else
								cumsum.erosion[slope.mainStation] = vecXdata[slope.mainStation].ErosionMass;
						} else {
							// Potential erosion at main station only
							if (cumsum.erosion[slope.mainStation] < -Constants::eps)
								cumsum.erosion[slope.mainStation] -= surfFluxes.mass[SurfaceFluxes::MS_WIND];
							else if (!(cumsum.erosion[slope.mainStation] > Constants::eps))
								cumsum.erosion[slope.mainStation] = -surfFluxes.mass[SurfaceFluxes::MS_WIND];
						}
					}
					const unsigned int i_hz = mn_ctrl.HzStep;
					if (mode == "OPERATIONAL") {
						if (!cumsum_mass) { // Cumulate flat field runoff in operational mode
							qr_Hdata.at(i_hz).runoff += surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];
							cumsum.runoff += surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];
						}
						/*
						 * Snow depth and mass corrections (deflate-inflate):
						 *   Monitor snow depth discrepancy assumed to be due to ...
						 *   ... wrong settling, which in turn is assumed to be due to a wrong estimation ...
						 *   of fresh snow mass because Michi spent many painful days calibrating the settling ...
						 *   and therefore it can't be wrong, dixunt Michi and Charles.
						 */
						const double cH = vecXdata[slope.mainStation].cH - vecXdata[slope.mainStation].Ground;
						const double mH = vecXdata[slope.mainStation].mH - vecXdata[slope.mainStation].Ground;
						// Look for missed erosion or not strong enough settling ...
						// ... and nastily deep "dips" caused by buggy data ...
						if (time_count_deltaHS > -Constants::eps2) {
							if ((mH + 0.01) < cH) {
								time_count_deltaHS += S_TO_D(sn_dt);
							} else {
								time_count_deltaHS = 0.;
							}
						}
						// ... or too strong settling
						if (time_count_deltaHS < Constants::eps2) {
							if ((mH - 0.01) > cH) {
								time_count_deltaHS -= S_TO_D(sn_dt);
							} else {
								time_count_deltaHS = 0.;
							}
						}
						// If the error persisted for at least one day => apply correction
						if (fabs(time_count_deltaHS) > (1. - 0.05 * M_TO_D(calculation_step_length))) {
							deflateInflate(Mdata, vecXdata[slope.mainStation],
							               qr_Hdata.at(i_hz).dhs_corr, qr_Hdata.at(i_hz).mass_corr);
							if (prn_check) {
								prn_msg(__FILE__, __LINE__, "msg+", Mdata.date,
								        "InflDefl (i_hz=%u): dhs=%f, dmass=%f, counter=%f",
								        i_hz, qr_Hdata.at(i_hz).dhs_corr, qr_Hdata.at(i_hz).mass_corr,
								        time_count_deltaHS);
							}
							time_count_deltaHS = 0.;
						}
					}
					if (mn_ctrl.HzDump) { // Save hazard data ...
						strncpy(qr_Hdata.at(i_hz).stat_abbrev, vecStationIDs[i_stn].c_str(), 15);
						if (mode == "OPERATIONAL") {
							qr_Hdata.at(i_hz).loc_for_snow = (int)vecStationIDs[i_stn][vecStationIDs[i_stn].length()-1];
							//TODO: WHAT SHOULD WE SET HERE? wstat_abk (not existing yet in DB) and wstao_nr, of course;-)
							qr_Hdata_ind.at(i_hz).loc_for_wind = -1;
						} else {
							qr_Hdata.at(i_hz).loc_for_snow = 2;
							qr_Hdata.at(i_hz).loc_for_wind = 1;
						}
						hazard.getHazardDataMainStation(qr_Hdata.at(i_hz), qr_Hdata_ind.at(i_hz),
						                                sn_Zdata, cumsum.drift, slope.mainStationDriftIndex,
						                                vecXdata[slope.mainStation], Mdata, surfFluxes);
						mn_ctrl.HzStep++;
						if (slope.mainStationDriftIndex)
							cumsum.drift = 0.;
						surfFluxes.hoar = 0.;
					}
					// New snow water equivalent (kg m-2), rain was dealt with in Watertransport.cc
					surfFluxes.mass[SurfaceFluxes::MS_HNW] += vecXdata[slope.mainStation].hn
					                                              * vecXdata[slope.mainStation].rho_hn;
					if (!avgsum_time_series) { // Sum up precipitations
						cumsum.rain += surfFluxes.mass[SurfaceFluxes::MS_RAIN];
						cumsum.snow += surfFluxes.mass[SurfaceFluxes::MS_HNW];
					}
				} else if (slope.luvDriftIndex) {
					const unsigned int i_hz = (mn_ctrl.HzStep > 0) ? mn_ctrl.HzStep-1 : 0;
					// Update drifting snow index (VI24),
					//   considering only snow eroded from the windward slope
					cumulate(cumsum.drift, surfFluxes.drift);
					if (mn_ctrl.HzDump) {
						// NOTE qr_Hdata was first saved at the end of the mainStation simulation, at which time the drift index could not be computed!
						hazard.getHazardDataSlope(qr_Hdata.at(i_hz), qr_Hdata_ind.at(i_hz),
						                          sn_Zdata.drift24, cumsum.drift, vecXdata[slope.sector], 
						                          slope.luvDriftIndex, slope.north, slope.south);
						cumsum.drift = 0.;
					}
					
					// Update erosion mass from windward virtual slope
					cumsum.erosion[slope.sector] += vecXdata[slope.sector].ErosionMass;
				}

				// TIME SERIES (*.met)
				if (tswrite && mn_ctrl.TsDump) {
					// Average fluxes
					if (avgsum_time_series) {
						averageFluxTimeSeries(mn_ctrl.nAvg, useCanopyModel, surfFluxes, vecXdata[slope.sector]);
					} else {
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] = cumsum.rain;
						surfFluxes.mass[SurfaceFluxes::MS_HNW] = cumsum.snow;
						// Add eroded snow from luv to precipitations on lee slope
						if (slope.sector == slope.lee && cumsum.erosion[slope.luv] > Constants::eps)
							surfFluxes.mass[SurfaceFluxes::MS_HNW] += cumsum.erosion[slope.luv] / vecXdata[slope.luv].cos_sl;
					}

					if (precip_rates) { // Precip rates in kg m-2 h-1
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
						surfFluxes.mass[SurfaceFluxes::MS_HNW] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
						if ((mode == "OPERATIONAL") && (!cumsum_mass)) {
							surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] = cumsum.runoff;
							surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
							cumsum.runoff = 0.;
						}
					}

					// Erosion mass rate in kg m-2 h-1
					surfFluxes.mass[SurfaceFluxes::MS_WIND] = cumsum.erosion[slope.sector];
					surfFluxes.mass[SurfaceFluxes::MS_WIND] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);

					// Dump
					const unsigned int i_hz = (mn_ctrl.HzStep > 0) ? mn_ctrl.HzStep - 1 : 0;
					unsigned int i_hz0 = (mn_ctrl.HzStep > 1) ? mn_ctrl.HzStep - 2 : 0;
					if (slope.mainStationDriftIndex)
						i_hz0 = i_hz;
					const double wind_trans24 = (slope.sector == slope.mainStation) ? qr_Hdata.at(i_hz0).wind_trans24 : qr_Hdata.at(i_hz).wind_trans24;
					snowpackio.writeTimeSeries(vecXdata[slope.sector], surfFluxes, Mdata,
					                           qr_Hdata.at(i_hz), wind_trans24);

					if (avgsum_time_series) {
						surfFluxes.reset(cumsum_mass);
						if (useCanopyModel) vecXdata[slope.sector].Cdata.reset(cumsum_mass);
					}
					surfFluxes.cRho_hn = Constants::undefined;
					surfFluxes.mRho_hn = Constants::undefined;
					// reset cumulative variables
					if (slope_sequence == slope.nSlopes-1) {
						cumsum.erosion.assign(cumsum.erosion.size(), 0.);
						cumsum.rain = cumsum.snow = 0.;
						mn_ctrl.nAvg = 0;
					}
				}

				// SNOW PROFILES ...
				// ... for visualization (*.pro), etc. (*.prf)
				if (profwrite && mn_ctrl.PrDump)
					snowpackio.writeProfile(current_date, vecXdata[slope.sector]);

				// ... backup Xdata (*.sno<JulianDate>)
				if (mn_ctrl.XdataDump) {
					std::stringstream ss;
					ss << vecStationIDs[i_stn];
					if (slope.sector != slope.mainStation) ss << slope.sector;
					snowpackio.writeSnowCover(current_date, vecXdata[slope.sector],
					                          vecSSdata[slope.sector], sn_Zdata, true);
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					        "Backup Xdata dumped for station %s [%.2f days, step %d]", ss.str().c_str(),
					        (current_date.getJulian()
					            - (vecSSdata[slope.mainStation].profileDate.getJulian() + 0.5/24)),
					        mn_ctrl.nStep);
				}

				// check mass balance if AVGSUM_TIME_SERIES is not set (screen output only)
				if (!avgsum_time_series) {
					const bool mass_balance = cfg.get("MASS_BALANCE", "SnowpackAdvanced");
					if (mass_balance) {
						if (massBalanceCheck(vecXdata[slope.sector], surfFluxes, tot_mass_in) == false)
							prn_msg(__FILE__, __LINE__, "msg+", current_date, "Mass error at end of time step!");
					}
				}
			} //end loop on slopes
			computed_one_timestep = true;
		} while ((dateEnd.getJulian() - current_date.getJulian()) > calculation_step_length/(2.*1440));
		//end loop on timesteps

		// If the simulation run for at least one time step,
		//   dump the PROFILEs (Xdata) for every station referred to as sector where sector 0 corresponds to the main station
		if (computed_one_timestep) {
			for (size_t sector=slope.mainStation; sector<slope.nSlopes; sector++) {
				if ((mode == "OPERATIONAL") && (sector == slope.mainStation)) {
					// Operational mode ONLY: dump snow depth discrepancy time counter
					vecSSdata[slope.mainStation].TimeCountDeltaHS = time_count_deltaHS;
				}
				snowpackio.writeSnowCover(current_date, vecXdata[sector], vecSSdata[sector], sn_Zdata);
				if (sector == slope.mainStation) {
					prn_msg(__FILE__, __LINE__, "msg", mio::Date(),
					        "Writing data to sno file(s) for %s (station %s) on %s",
					        vecSSdata[slope.mainStation].meta.getStationName().c_str(),
					        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str());
				}
			}
			// Dump time series to snowpack.ams_pmod@SDBx (hazard data)
			if (mn_ctrl.sdbDump) {
				mio::Timer sdbDump_timer;
				sdbDump_timer.reset();
				sdbDump_timer.start();
				if (snowpackio.writeHazardData(vecStationIDs[i_stn], qr_Hdata, qr_Hdata_ind, mn_ctrl.HzStep)) {
					sdbDump_timer.stop();
					prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
					        "Finished writing Hdata to SDB for station %s on %s (%lf s)",
					        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str(), sdbDump_timer.getElapsed());
				}
			}
		}
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Total time to read meteo data : %lf s",
		        meteoRead_timer.getElapsed());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Runtime for station %s: %lf s",
		        vecStationIDs[i_stn].c_str(), run_timer.getElapsed());
	}

	time_t nowEND=time(NULL);
	cout << endl;
	cout << "[i] []                 STARTED  running SLF " << mode << " Snowpack Model on " << ctime(&nowSRT);
	if (mode == "OPERATIONAL"){
		cout << "                       ===========================================================================" << endl;
	} else {
		cout << "                       ========================================================================" << endl;
	}
	cout << "                       FINISHED running SLF " << mode << " Snowpack Model on " << ctime(&nowEND) << endl;
}

int main(int argc, char *argv[]) {
	try {
		real_main(argc, argv);
	} catch (const std::exception &e) {
		std::cerr << e.what() << endl;
		exit(1);
	}

	return 0;
}
