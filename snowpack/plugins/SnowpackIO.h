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

#ifndef SNOWPACKIO_H
#define SNOWPACKIO_H

#include <meteoio/MeteoIO.h>

#include <snowpack/DataClasses.h>
#include <snowpack/plugins/SnowpackIOInterface.h>
#include <snowpack/plugins/SmetIO.h>
#include <snowpack/plugins/AsciiIO.h>
#ifdef CAAMLIO
	#include <snowpack/plugins/CaaMLIO.h>
#endif
#ifdef IMISDBIO
	#include <snowpack/plugins/ImisDBIO.h>
#endif

/**
 * @page snowpackio Snowpack data formats
 * 
 * @section snowpack_inputs_outputs Snowpack inputs versus outputs
 * %Snowpack needs several kind of information to be provided for a simulation and then writes out
 * several kind of information. Some formats can be used for both the inputs and the outputs while some
 * others are restricted to either one or the other (simply because %Snowpack does not read or write 
 * out the parameters contained in the said format).
 * 
 * @subsection Snowpack_inputs Snowpack required inputs
 *  Several kind of information need to be given to %Snowpack for a simulation:
 * -# the description of the place where the snow pack has to be simulated: latitutde, longitude, elevation, slope, ...
 * -# the time series of the various meteorological parameters
 * -# the intial state of the various soil and snow layers
 *
 * Very often, 1) and 2) are provided together. But this depends ultimately on the file format that is used ot provide such data (SMET, INP, etc). These two points are
 * handled by <a href="https://models.slf.ch/p/meteoio">MeteoIO</a>, so please check its documentation, in the 
 * <i>"Available plugins and usage"</i> section for the relevant formats.
 * 
 * @subsection Snowpack_outputs Snowpack outputs
 * %Snowpack creates various output files:
 * - the current state of its soil and snow layers in <i>".sno"</i> files;
 * - the current state of its hazard relevant data in <i>".haz"</i> files;
 * - a time serie of snow profiles;
 * - a time serie of the meteorological data and fluxes as used in the model.
 * 
 * Depending on the chosen output format, 1) and 2) might be provided as one file or two files.
 *
 * @section available_single_profile_plugins Single snow profiles
 * The %Snowpack specific data are supported directly in %Snowpack and the formats listed in the table below 
 * are available, both for input and output of snow profiles with the <b>"SNOW"</b> keyword. 
 * Please read the documentation for each plugin in order to know the plugin-specific keywords!
 * <center><table border="1">
 * <tr><th>Key</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage snoold_format "SNOOLD"</td><td>legacy %Snowpack profile (including the hazard data)</td><td></td></tr>
 * <tr><td>\subpage smet "SMET"</td><td>SMET based profile (including the hazard data), recommended</td><td></td></tr>
 * <tr><td>\subpage caaml "CAAML"</td><td>CAAML profile</td><td><A HREF="http://xmlsoft.org/">libxml</A></td></tr>
 * </table></center>
 * 
 * @section available_profile_ts_plugins Snow profiles time series
 * The %Snowpack specific data are supported directly in %Snowpack and the formats listed in the table below 
 * are available for output of snow profiles time series with the <b>"PROFILE_FORMAT"</b> keyword. 
 * Please read the documentation for each plugin in order to know the plugin-specific keywords!
 * <center><table border="1">
 * <tr><th>Key</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage pro_format "PRO"</td><td>legacy %Snowpack profile time series</td><td></td></tr>
 * <tr><td>\subpage prf_format "PRF"</td><td>easier to parse profile time series</td><td></td></tr>
 * <tr><td>\subpage profile_imis "IMIS"</td><td>write profile time series to the IMIS database</td><td><A HREF="http://docs.oracle.com/cd/B12037_01/appdev.101/b10778/introduction.htm">Oracle's OCCI library</A></td></tr>
 * </table></center>
 * 
 * @section available_met_ts Fluxes time series
 * %Snowpack computes various meteorological parameters as well as fluxes and can write them out as time series.
 * Currently, only the \subpage met_format "MET format" is supported.
 *
 */
class SnowpackIO : public SnowpackIOInterface {

	public:
		SnowpackIO(const SnowpackConfig& cfg);
		SnowpackIO(const SnowpackIO& source);
		~SnowpackIO();

		virtual bool snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const;

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
                                    const ZwischenData& Zdata, const bool& forbackup=false);

		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);

		virtual void writeProfile(const mio::Date& date, const SnowStation& Xdata);

		virtual bool writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                     const std::vector<ProcessInd>& Hdata_ind, const size_t& num);

		SnowpackIO& operator=(const SnowpackIO& source);

	private:
#ifdef IMISDBIO
		ImisDBIO *imisdbio;
#endif
#ifdef CAAMLIO
		CaaMLIO *caamlio;
#endif
		SmetIO *smetio;
		AsciiIO *asciiio;
		bool input_snow_as_smet, output_snow_as_smet;
		bool input_snow_as_caaml, output_snow_as_caaml;
		bool input_snow_as_ascii, output_snow_as_ascii;
		bool output_prf_as_ascii, output_prf_as_caaml, output_prf_as_imis;
		bool output_ts_as_ascii, output_haz_as_imis;
};

#endif //End of SnowpackIO.h
