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

#include <snowpack/plugins/SnowpackIO.h>

#include <snowpack/plugins/SmetIO.h>
#include <snowpack/plugins/AsciiIO.h>

#cmakedefine PLUGIN_IMISIO
#cmakedefine PLUGIN_CAAMLIO

#ifdef PLUGIN_CAAMLIO
	#include <snowpack/plugins/CaaMLIO.h>
#endif
#ifdef PLUGIN_IMISIO
	#include <snowpack/plugins/ImisDBIO.h>
#endif

using namespace std;
using namespace mio;

SnowpackIO::SnowpackIO(const SnowpackConfig& cfg):
	imisdbio(NULL), caamlio(NULL), smetio(NULL), asciiio(NULL),
	input_snow_as_smet(false), output_snow_as_smet(false),
	input_snow_as_caaml(false), output_snow_as_caaml(false),
	input_snow_as_ascii(false), output_snow_as_ascii(false),
	output_prf_as_ascii(false), output_prf_as_caaml(false), output_prf_as_imis(false),
	output_ts_as_ascii(false), output_ts_as_smet(false), output_haz_as_imis(false)

{
	//Format of initial snow profile:
	//TODO: document Input::SNOW = SMET, CAAML, or SNOOLD
	const string in_snow = cfg.get("SNOW", "Input", IOUtils::nothrow);
	if (in_snow == "SNOOLD") {
		input_snow_as_ascii = true;
	} else if (in_snow == "CAAML") {
		input_snow_as_caaml = true;
	} else if (in_snow == "SMET") {
		input_snow_as_smet = true;
	} else
		throw InvalidArgumentException("Invalid input snow profile format '"+in_snow+"'. Please choose from SMET, CAAML, SNOOLD", AT);

	//Format of transitional and final snow profile(s):
	//TODO: document ouput::SNOW = SMET, CAAML, or SNOOLD
	const string out_snow = cfg.get("SNOW", "Output", IOUtils::nothrow);
	if (out_snow == "SNOOLD") {
		output_snow_as_ascii = true;
	} else if (out_snow == "CAAML") {
		output_snow_as_caaml = true;
	} else if (out_snow == "SMET") {
		output_snow_as_smet = true;
	} else
		throw InvalidArgumentException("Invalid output snow profile format '"+out_snow+"'. Please choose from SMET, CAAML, SNOOLD", AT);

	/* Profiles may be dumped in up to 3 formats specified by the key PROF_FORMAT in [Output].
	 * Note that the keys AGGREGATE_PRO and AGGREGATE_PRF will allow to aggregate model layers to a smaller number.
	 * PRO   : Full profiles in ASCII-format, including soil elements if available, for visualization with both SnopViz and SN_GUI
	 * PRF   : Snow profiles in tabular ASCII-format
	 * IMIS  : aggregated snow profiles for upload to the SLF database sdbo
	 */
	std::vector<string> vecProfileFmt = cfg.get("PROF_FORMAT", "Output", IOUtils::nothrow);
	if (vecProfileFmt.size() > 3) {
		throw InvalidArgumentException("The key PROF_FORMAT in [Output] can take three values at most", AT);
	} else {
		for (size_t ii=0; ii<vecProfileFmt.size(); ii++) {
			if (vecProfileFmt[ii] == "PRO") {
				output_prf_as_ascii = true;
			} else if (vecProfileFmt[ii] == "PRF") {
				output_prf_as_ascii  = true;
			} else if (vecProfileFmt[ii] == "IMIS") {
				output_prf_as_imis  = true;
			} else {
				throw InvalidArgumentException("The key PROF_FORMAT in [Output] takes only PRO, PRF or IMIS as value", AT);
			}
		}
	}

	//Format of meteo time series:
	const bool ts_out = cfg.get("TS_WRITE", "Output");
	if (ts_out==true) {
		const std::string ts_format = cfg.get("TS_FORMAT", "Output", IOUtils::nothrow);
		if (ts_format=="SMET")
			output_ts_as_smet = true;
		else if (ts_format=="MET")
			output_ts_as_ascii = true;
		else
			throw InvalidArgumentException("The key TS_FORMAT in [Output] takes only SMET or MET as value", AT);
	}

	//set the "plugins" pointers
	RunInfo run_info;
	if (input_snow_as_smet || output_snow_as_smet || output_ts_as_smet) smetio = new SmetIO(cfg, run_info);
	if (input_snow_as_ascii || output_snow_as_ascii || output_prf_as_ascii || output_ts_as_ascii) asciiio = new AsciiIO(cfg, run_info);
#ifdef PLUGIN_CAAMLIO
	if (input_snow_as_caaml || output_snow_as_caaml) caamlio = new CaaMLIO(cfg, run_info);
#endif
#ifdef PLUGIN_IMISIO
	output_haz_as_imis = output_prf_as_imis;
	if (output_prf_as_imis || output_haz_as_imis) imisdbio = new ImisDBIO(cfg, run_info);
#endif
}

SnowpackIO::SnowpackIO(const SnowpackIO& source) :
	imisdbio(source.imisdbio), caamlio(source.caamlio), smetio(source.smetio), asciiio(source.asciiio),
	input_snow_as_smet(source.input_snow_as_smet), output_snow_as_smet(source.input_snow_as_smet),
	input_snow_as_caaml(source.input_snow_as_caaml), output_snow_as_caaml(source.output_snow_as_caaml),
	input_snow_as_ascii(source.input_snow_as_ascii), output_snow_as_ascii(source.output_snow_as_ascii),
	output_prf_as_ascii(source.output_prf_as_ascii), output_prf_as_caaml(source.output_prf_as_caaml), output_prf_as_imis(source.output_prf_as_imis),
	output_ts_as_ascii(source.output_ts_as_ascii), output_ts_as_smet(source.output_ts_as_smet), output_haz_as_imis(source.output_haz_as_imis)
{}

SnowpackIO::~SnowpackIO()
{
	if (smetio != NULL) delete smetio;
	if (asciiio != NULL) delete asciiio;
	if (caamlio != NULL) delete caamlio;
	if (imisdbio != NULL) delete imisdbio;
}

bool SnowpackIO::snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const
{
	if (input_snow_as_ascii) {
		return asciiio->snowCoverExists(i_snowfile, stationID);
#ifdef PLUGIN_CAAMLIO
	} else if (input_snow_as_caaml){
		return caamlio->snowCoverExists(i_snowfile, stationID);
#endif
	} else {
		return smetio->snowCoverExists(i_snowfile, stationID);
	}
}

void SnowpackIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                               SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata)
{
	if (input_snow_as_ascii) {
		asciiio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
#ifdef PLUGIN_CAAMLIO
	} else if (input_snow_as_caaml) {
		caamlio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
#endif
	} else {
		smetio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	}
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
                                const ZwischenData& Zdata, const bool& forbackup)
{
	if (output_snow_as_ascii) {
		asciiio->writeSnowCover(date, Xdata, Zdata, forbackup);
#ifdef PLUGIN_CAAMLIO
	} else if (output_snow_as_caaml) {
		caamlio->writeSnowCover(date, Xdata, Zdata, forbackup);
#endif
	} else {
		smetio->writeSnowCover(date, Xdata, Zdata, forbackup);
	}
}

void SnowpackIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
                                 const ProcessDat& Hdata, const double wind_trans24)
{
	if (output_ts_as_ascii)
		asciiio->writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
	else if (output_ts_as_smet)
		smetio->writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
}

void SnowpackIO::writeProfile(const mio::Date& date, const SnowStation& Xdata)
{
	if (output_prf_as_ascii)
		asciiio->writeProfile(date, Xdata);

#ifdef PLUGIN_CAAMLIO
	if (output_prf_as_caaml)
		caamlio->writeProfile(date, Xdata);
#endif

#ifdef PLUGIN_IMISIO
	if (output_prf_as_imis)
		imisdbio->writeProfile(date, Xdata);
#endif
}

#ifdef PLUGIN_IMISIO
bool SnowpackIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const size_t& num)
{
	if(output_haz_as_imis)
		return imisdbio->writeHazardData(stationID, Hdata, Hdata_ind, num);
	return false;
}
#else
bool SnowpackIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                                 const std::vector<ProcessInd>& /*Hdata_ind*/, const size_t& /*num*/)
{
	return false;
}
#endif

SnowpackIO& SnowpackIO::operator=(const SnowpackIO& source)
{
	if(this != &source) {
		imisdbio = source.imisdbio;
		caamlio = source.caamlio;
		asciiio = source.asciiio;
		smetio = source.smetio;
		output_prf_as_ascii = source.output_prf_as_ascii;
		output_prf_as_caaml = source.output_prf_as_caaml;
		output_prf_as_imis = source.output_prf_as_imis;
		output_snow_as_caaml = source.output_snow_as_caaml;
		output_snow_as_smet = source.output_snow_as_smet;
		input_snow_as_caaml = source.input_snow_as_caaml;
		input_snow_as_smet = source.input_snow_as_smet;
		output_ts_as_ascii = source.output_ts_as_ascii;
		output_haz_as_imis = source.output_haz_as_imis;
	}
	return *this;
}
