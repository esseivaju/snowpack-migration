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

#include <snowpack/SnowpackIO.h>

using namespace std;
using namespace mio;

#ifdef IMISDBIO
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : imisdbio(i_cfg),
#else
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) :
#endif
              asciiio(i_cfg), smetio(i_cfg),
              cfg(i_cfg), outputprofile_as_ascii(false), outputprofile_as_imis(false),
              output_snow_as_smet(false), input_snow_as_smet(false)
{
	//The profiles may be dumped either in ASCII format or in another ASCII format for upload to the DB
	//The user can switch the desired mode on by specifying "ASCII" or "IMIS" or both in the io.ini
	vector<string> vecProfileOutput = cfg.get("PROFILE", "Output", IOUtils::nothrow);
	if (vecProfileOutput.empty()) {
		outputprofile_as_ascii = true;
		outputprofile_as_imis  = false;
	} else if (vecProfileOutput.size() > 2) {
		throw InvalidArgumentException("The key PROFILE in section OUTPUT can have two values at most", AT);
	} else {
		for (unsigned int ii=0; ii<vecProfileOutput.size(); ii++){
			if (vecProfileOutput[ii] == "ASCII"){
				outputprofile_as_ascii = true;
			} else if (vecProfileOutput[ii] == "IMIS") {
				outputprofile_as_imis  = true;
			} else {
				throw InvalidArgumentException("Key PROFILE / section OUTPUT: only values ASCII or IMIS expected", AT);
			}
		}
	}

	//Format of initial snow profile:
	const string in_snow = cfg.get("SNOW", "Input", IOUtils::nothrow);
	if (in_snow == "SMET"){ //TODO: document Input::SNOW = SMET or SNOOLD
		input_snow_as_smet = true;
	}
	//Format of transitional and final snow profile(s):
	const string out_snow = cfg.get("SNOW", "Output", IOUtils::nothrow);
	if (out_snow == "SMET"){ //TODO: document ouput::SNOW = SMET or SNOOLD
		output_snow_as_smet = true;
	}
}

bool SnowpackIO::snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const
{
	if (input_snow_as_smet){
		return smetio.snowCoverExists(i_snowfile, stationID);
	} else {
		return asciiio.snowCoverExists(i_snowfile, stationID);
	}
}

void SnowpackIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                               SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata)
{
	if (input_snow_as_smet){
		smetio.readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	} else {
		asciiio.readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	}
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
                                const ZwischenData& Zdata, const bool& forbackup)
{
	if (output_snow_as_smet){
		smetio.writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
	} else {
		asciiio.writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
	}
}

void SnowpackIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
                                 const ProcessDat& Hdata, const double wind_trans24)
{
	asciiio.writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
}

void SnowpackIO::writeProfile(const mio::Date& date, SnowStation& Xdata, const ProcessDat& Hdata)
{
	if (outputprofile_as_ascii)
		asciiio.writeProfile(date, Xdata, Hdata);

	if (outputprofile_as_imis){
#ifdef IMISDBIO
		imisdbio.writeProfile(date, Xdata, Hdata);
#endif
	}
}

#ifdef IMISDBIO
bool SnowpackIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const int& num)
{
	if (imisdbio.writeHazardData(stationID, Hdata, Hdata_ind, num)){
		return true;
	} else {
		return false;
	}
}
#else
bool SnowpackIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                                 const std::vector<ProcessInd>& /*Hdata_ind*/, const int& /*num*/)
{
	return false;
}
#endif
