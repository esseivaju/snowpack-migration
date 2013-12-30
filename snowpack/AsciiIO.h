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

#ifndef __ASCIIIO_H__
#define __ASCIIIO_H__

#include <snowpack/Constants.h>
#include <snowpack/SnowpackIOInterface.h>
#include <snowpack/Hazard.h>
#include <meteoio/MeteoIO.h>
#include <snowpack/Canopy.h>

class AsciiIO : public SnowpackIOInterface {

	public:
		AsciiIO(const SnowpackConfig& i_cfg, const RunInfo& run_info);

		virtual bool snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const;

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
		                            const ZwischenData& Zdata, const bool& forbackup=false);

		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);

		virtual void writeProfile(const mio::Date& date, const SnowStation& Xdata);

		virtual bool writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
		                             const std::vector<ProcessInd>& Hdata_ind, const int& num);

	private:
		typedef enum {
			PRO,  ///< Profile visualization file
			FULL, ///< Full tabular ascii profile, 1 element == 1 layer
			AGGR  ///< Aggregated tabular ascii profile, N elements == 1 layer
		} PRF_TYPE;

		bool appendFile(const std::string& filename, const mio::Date& startdate, const std::string& ftype);
		bool parseMetFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp);
		bool parseProFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp);

		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true) const;

		bool checkHeader(const SnowStation& Xdata, const std::string& filename, const std::string& ext, const std::string& signature) const;

		void writeProfilePro(const mio::Date& date, const SnowStation& Xdata);
		void writeProfileProAddDefault(const SnowStation& Xdata, FILE *fout);
		void writeProfileProAddCalibration(const SnowStation& Xdata, FILE *fout);

		void writeProfilePrf(const mio::Date& date, const SnowStation& Xdata, const bool& aggregate);

		size_t writeTemperatures(FILE *fout, const double& z_vert, const double& T,
		                         const size_t& ii, const SnowStation& Xdata);

		double compPerpPosition(const double& z_vert, const double& hs_ref,
		                        const double& ground, const double& cos_sl);
		double checkMeasuredTemperature(const double& T, const double& z, const double& mH);

		size_t findTaggedElement(const size_t& tag, const SnowStation& Xdata);
		size_t writeHeightTemperatureTag(FILE *fout, const size_t& tag,
		                                 const CurrentMeteo& Mdata, const SnowStation& Xdata);

		void setNumberSensors(const CurrentMeteo& Mdata);
		void writeTimeSeriesAddDefault(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                       const CurrentMeteo& Mdata, const double crust,
                                       const double dhs_corr, const double mass_corr,
                                       const size_t nCalcSteps, FILE *fout);
		void writeTimeSeriesAddAntarctica(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                          const CurrentMeteo& Mdata, const double crust,
                                          const double dhs_corr, const double mass_corr,
                                          const size_t nCalcSteps, FILE *fout);
		void writeTimeSeriesAddCalibration(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                           const CurrentMeteo& Mdata, const double crust,
                                           const double dhs_corr, const double mass_corr,
                                           const size_t nCalcSteps, FILE *fout);

		void readTags(const std::string& filename, const CurrentMeteo&  Mdata, TaggingData& TAGdata);

		std::set<std::string> setAppendableFiles;
		std::string variant, experiment, sw_mode;
		std::string inpath, snowfile, i_snopath, outpath, o_snopath;
		const RunInfo info;

		std::vector<std::string> vecProfileFmt;
		bool aggregate_prf;

		//Monitored temperature sensors
		std::vector<double> fixedPositions;
		size_t numberMeasTemperatures, maxNumberMeasTemperatures;
		size_t numberTags, numberFixedSensors, totNumberSensors;

		double time_zone; // input data time zone
		double calculation_step_length, hazard_steps_between, ts_days_between;
		double min_depth_subsurf, hoar_density_surf, hoar_min_size_surf;
		bool avgsum_time_series, useCanopyModel, useSoilLayers, research_mode, perp_to_slope;
		bool out_heat, out_lw, out_sw, out_meteo, out_haz, out_mass, out_t, out_load, out_stab, out_canopy;
		bool r_in_n;

		static const bool t_srf, t_gnd;
};

#endif //End of AsciiIO.h
