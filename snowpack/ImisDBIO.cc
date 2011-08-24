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

#include <snowpack/ImisDBIO.h>

using namespace std;
using namespace mio;
using namespace oracle;
using namespace oracle::occi;

const double ImisDBIO::time_zone = 1.; //All IMIS data is in gmt+1

double ImisDBIO::hoar_density_surf = 0.0;
double ImisDBIO::hoar_min_size_surf = 0.0;

const string ImisDBIO::sqlDeleteHdata = "DELETE FROM snowpack.ams_pmod WHERE stat_abk=:1 and stao_nr=:2 and datum>=:3 and datum<=:4";

const string ImisDBIO::sqlInsertHdata = "INSERT INTO snowpack.ams_pmod(datum,stat_abk,stao_nr,dewpt_def,hoar_ind6,hoar_ind24,wind_trans,hns3,hns6,hns12,hns24,hns72,hns72_24,wc3,wc6,wc12,wc24,wc72,hoar_size,wind_trans24,stab_class1,stab_class2,stab_index1,stab_height1,stab_index2,stab_height2,stab_index3,stab_height3,stab_index4,stab_height4,stab_index5,stab_height5,ch,crust,en_bal,sw_net,t_top1,t_top2,snowpack_version,calc_date,swe,tot_lwc,runoff) values (:1,:2,:3,:4,:5,:6,:7,:8,:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:25,:26,:27,:28,:29,:30,:31,:32,:33,:34,:35,:36,:37,:38,:39,:40,:41,:42,:43)";

string ImisDBIO::oracleDB = "";
string ImisDBIO::oracleUser = "";
string ImisDBIO::oraclePassword = "";

const std::string ImisDBIO::profile_filename = "loaddata/pmodpro.dat";

ImisDBIO::ImisDBIO(const mio::Config& cfg)
{
	cfg.getValue("DBNAME", "Output", oracleDB, Config::nothrow);
	cfg.getValue("DBUSER", "Output", oracleUser, Config::nothrow);
	cfg.getValue("DBPASS", "Output", oraclePassword, Config::nothrow);

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoar_density_surf);

	//Minimum size to show surface hoar on surface (mm)
	cfg.getValue("HOAR_MIN_SIZE_SURF", "SnowpackAdvanced", hoar_min_size_surf);
}

void ImisDBIO::readSnowCover(const std::string& /*i_snowfile*/, const std::string& /*stationID*/,
                             SN_SNOWSOIL_DATA& /*SSdata*/, SN_ZWISCHEN_DATA& /*Zdata*/)
{
	throw IOException("Nothing implemented here!", AT);
}

//void ImisDBIO::writeSnowCover(const mio::Date& /*date*/, const std::string& /*station*/, const SnowStation& /*Xdata*/,
void ImisDBIO::writeSnowCover(const mio::Date& /*date*/, const SnowStation& /*Xdata*/,
                              const SN_SNOWSOIL_DATA& /*SSdata*/, const SN_ZWISCHEN_DATA& /*Zdata*/, const bool& /*forbackup*/)
{
	throw IOException("Nothing implemented here!", AT);
}

void ImisDBIO::writeTimeSeries(const SnowStation& /*Xdata*/, const SurfaceFluxes& /*Sdata*/, const CurrentMeteo& /*Mdata*/,
                               const ProcessDat& /*Hdata*/, const double /*wind_trans24*/)
{
	throw IOException("Nothing implemented here!", AT);
}

/**
 * @brief Dump snow profile to ASCII file for subsequent upload to SDBO
 */
void ImisDBIO::writeProfile(const mio::Date& date, SnowStation& Xdata, const ProcessDat& Hdata)
{
	unsigned int nE = Xdata.getNumberOfElements();
	if ((Xdata.sector != 0) || (nE == 0)) {
		return;
	}

	FILE *PFile=NULL;
	unsigned int e, ll;
	const vector<NodeData>& NDS = Xdata.Ndata;
	const vector<ElementData>& EMS = Xdata.Edata;

	vector<SnowProfileLayer> Pdata;
	Pdata.resize(nE);

	// Generate the profile data from the element data (1 layer = 1 element)
	int snowloc = 0;
	string mystation = Xdata.meta.getStationID();
	if (isdigit(mystation[mystation.length()-1])) {
		snowloc = mystation[mystation.length()-1] - '0';
		if (mystation.length() > 2)
			mystation = mystation.substr(0, mystation.length()-1);
	}

	for(e=0, ll=0; e<nE; e++) {
		Pdata[ll].profileDate = date;
		Pdata[ll].stationname = mystation;
		Pdata[ll].loc_for_snow = snowloc;
		Pdata[ll].loc_for_wind = 1;

		//write version and date
		Pdata[ll].sn_version = Hdata.sn_version;
		Pdata[ll].sn_computation_date = Hdata.sn_computation_date;
		Pdata[ll].sn_jul_computation_date=Hdata.sn_jul_computation_date;
		Pdata[ll].height = M_TO_CM(NDS[e+1].z + NDS[e+1].u);
		Pdata[ll].layerDate = EMS[e].depositionDate;
		Pdata[ll].rho = EMS[e].Rho;
		Pdata[ll].T = K_TO_C(NDS[e+1].T);
		Pdata[ll].gradT = EMS[e].gradT;
		Pdata[ll].strain_rate = fabs(EMS[e].EDot);
		Pdata[ll].theta_w = EMS[e].theta[WATER] * 100.;
		Pdata[ll].theta_i = EMS[e].theta[ICE] * 100.;
		Pdata[ll].dendricity = EMS[e].dd;
		Pdata[ll].sphericity = EMS[e].sp;
		Pdata[ll].coordin_num = EMS[e].N3;
		Pdata[ll].grain_size = 2. * EMS[e].rg;
		Pdata[ll].bond_size = 2. * EMS[e].rb;
		Pdata[ll].marker = EMS[e].mk%100;
		Pdata[ll].hard = EMS[e].hard;
		ll++;
	}

	unsigned int nL = Aggregate::aggregate(Pdata);

	if (!(PFile = fopen(profile_filename.c_str(), "a"))) {
		prn_msg(__FILE__, __LINE__, "err", date, "Cannot open Profile file: %s", profile_filename.c_str());
		throw FileAccessException(profile_filename, AT);
	}

	for(ll=0; ll<nL; ll++) {
		//HACK: these legacy offset should be removed.
		//This means specify a different import date format for the database and remove the offset here
		const double profile_date = Pdata[ll].profileDate.getJulianDate() - 2415021. + 0.5; //HACK
		const double layer_date = Pdata[ll].layerDate.getJulianDate() - 2415021. + 0.5; //HACK

		fprintf(PFile,"%.5f,%s,%d,%.2f,", profile_date, Pdata[ll].stationname.c_str(),
		        Pdata[ll].loc_for_snow, Pdata[ll].height);
		fprintf(PFile,"%.5f,%.0f,%.1f,%.0f,%.4e,%.0f,%.0f,%.2f,%.2f,%.1f,%.1f,%.2f,%d\n",
		        layer_date, Pdata[ll].rho, Pdata[ll].T, Pdata[ll].gradT, Pdata[ll].strain_rate,
		        Pdata[ll].theta_w, Pdata[ll].theta_i, Pdata[ll].dendricity, Pdata[ll].sphericity,
		        Pdata[ll].coordin_num, Pdata[ll].grain_size, Pdata[ll].bond_size, Pdata[ll].type);
	}

	if (NDS[nE].hoar > MM_TO_M(hoar_min_size_surf) * hoar_density_surf) {
		ll=nL-1;
		//HACK: these legacy offset should be removed.
		//This means specify a different import date format for the database and remove the offset here
		const double profile_date = Pdata[ll].profileDate.getJulianDate() - 2415021. + 0.5; //HACK
		const double layer_date = Pdata[ll].layerDate.getJulianDate() - 2415021. + 0.5; //HACK
		double gsz_SH = NDS[nE].hoar / hoar_density_surf;
		const double Tss = Pdata[ll].T + (Pdata[ll].gradT * gsz_SH);

		fprintf(PFile,"%.5f,%s,%d,%.2f,", profile_date, Pdata[ll].stationname.c_str(),
		        Pdata[ll].loc_for_snow, Pdata[ll].height + M_TO_CM(gsz_SH));
		fprintf(PFile,"%.5f,%.0f,%.1f,%.0f,%.4e,%.0f,%.0f,%.2f,%.2f,%.1f,%.1f,%.2f,%d\n",
		        layer_date, hoar_density_surf, Tss, Pdata[ll].gradT, 0.,
		        0., hoar_density_surf/Constants::density_ice, 0., 0., 2., M_TO_MM(gsz_SH), 0.6667*M_TO_MM(gsz_SH), 660);
	}

	fclose(PFile);
}

bool ImisDBIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                               const std::vector<ProcessInd>& Hdata_ind, const int& num)
{
	//HACK: num is incremented after each new data is added. It is therefore the index of the next element to write
	if ((num <= 0) || (num > (int)Hdata.size())){
		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "No hazard data either deleted from or inserted into %s: %d steps while Hdata.size=%d", oracleDB.c_str(), num, Hdata.size());
		return false; //nothing to do
	}

	if ((oracleDB == "") || (oraclePassword == "") || (oracleUser == "")){
		//throw IOException("You must set the output database, username and password", AT);
		if (num >= (int)Hdata.size()){
			prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "No data written to %s!", oracleDB.c_str());
			prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "You must set all of output database, username and password first");
		}
		return false; //nothing to do
	}

	string stationName="", stationNumber="";
	parseStationName(stationID, stationName, stationNumber);

	Environment *env = NULL;

	try {
		env = Environment::createEnvironment();// static OCCI function
		Connection *conn = NULL;

		conn = env->createConnection(oracleUser, oraclePassword, oracleDB);

		deleteHdata(stationName, stationNumber, Hdata[0].date, Hdata[num-1].date, env, conn);
		insertHdata(stationName, stationNumber, Hdata, Hdata_ind, num, env, conn);

		env->terminateConnection(conn);
		Environment::terminateEnvironment(env); // static OCCI function

	} catch (const exception& e){
		Environment::terminateEnvironment(env); // static OCCI function
		prn_msg(__FILE__, __LINE__, "err", mio::Date(), ":");
		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "while writing hazard data for %s to %s,",
		        stationID.c_str(), oracleDB.c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "from %s to %s",
		        Hdata[0].date.toString(mio::Date::ISO).c_str(), Hdata[num-1].date.toString(mio::Date::ISO).c_str());
		throw IOException("Oracle Error: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
	return true;
}

/**
 * @brief This function breaks up the station name into two components (a string and a number e.g. KLO2 -> "KLO","2")
 * @param stationName The full name of the station (e.g. "KLO2")
 * @param stName      The string part of the name  (e.g. "KLO")
 * @param stNumber    The integer part of the name (e.g. "2")
 */
void ImisDBIO::parseStationName(const std::string& stationName, std::string& stName, std::string& stNumber)
{
	stName    = stationName.substr(0, stationName.length()-1); //The station name: e.g. KLO
	stNumber  = stationName.substr(stationName.length()-1, 1); //The station number: e.g. 2
	if(!std::isdigit(stNumber[0])) {
		//the station is one of these non-imis stations that don't contain a number...
		stName = stationName;
		stNumber = "0";
	}
}

void ImisDBIO::deleteHdata(const std::string& stationName, const std::string& stationNumber,
                           const mio::Date& dateStart, const mio::Date& dateEnd,
                           oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	vector< vector<string> > vecResult;
	vector<int> datestart = vector<int>(5);
	vector<int> dateend   = vector<int>(5);

	//IMIS is in TIME_ZONE=+1, so moving back to this time_zone
	mio::Date dateS(dateStart), dateE(dateEnd);
	dateS.setTimeZone(time_zone);
	dateE.setTimeZone(time_zone);
	dateS.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	dateE.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);

	//Oracle can't deal with an integer for the hour of 24, hence the following workaround
	if (datestart[3] == 24){
		mio::Date tmpDate = dateS + 3.0/(60*60*24); //add three seconds to omit 24 for 00
		tmpDate.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	}
	if (dateend[3] == 24){
		mio::Date tmpDate = dateEnd + 3.0/(60*60*24); //add three seconds to omit 24 for 00
		tmpDate.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
	}

	Statement *stmt  = conn->createStatement(sqlDeleteHdata);

	// construct the oracle specific Date object: year, month, day, hour, minutes
	occi::Date begindate(env, datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	occi::Date enddate(env, dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
	stmt->setString(1, stationName);   // set 1st variable's value (station name)
	stmt->setString(2, stationNumber); // set 2nd variable's value (station number)
	stmt->setDate(3, begindate);       // set 4rd variable's value (begin date)
	stmt->setDate(4, enddate);         // set 5th variable's value (enddate)

	unsigned int rows_deleted = stmt->executeUpdate();
	conn->terminateStatement(stmt);

	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Deleted %d rows in %s!", rows_deleted, oracleDB.c_str());
}

void ImisDBIO::insertHdata(const std::string& stationName, const std::string& stationNumber,
                           const std::vector<ProcessDat>& Hdata, const std::vector<ProcessInd>& Hdata_ind,
                           const int& num, oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	vector<int> sndate = vector<int>(5);
	unsigned int rows_inserted = 0;
	double sn_version;
	IOUtils::convertString(sn_version, Hdata[0].sn_version);

	mio::Date dateSn( Hdata[0].sn_jul_computation_date, time_zone );
	dateSn.getDate(sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);
	if (sndate[3] == 24){
		mio::Date tmpDate = dateSn + 3.0/(60*60*24); //add three seconds to omit 24 for 00
		tmpDate.getDate(sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);
	}

	for (int i = 0; i<num; i++){
		if (Hdata[i].date == mio::Date()) break; //catch the case that not all Hdata has been set properly

		vector<int> hzdate = vector<int>(5);
		mio::Date dateH(Hdata[i].date);
		dateH.setTimeZone(time_zone);
		dateH.getDate(hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);

		//Oracle can't deal with an integer for the hour of 24, hence the following workaround
		if (hzdate[3] == 24){
			mio::Date tmpDate = dateH + 3.0/(60*60*24); //add three seconds to omit 24 for 00
			tmpDate.getDate(hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);
		}

		Statement *stmt  = conn->createStatement(sqlInsertHdata);

		// construct the oracle specific Date object: year, month, day, hour, minutes
		occi::Date hazarddate(env, hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);
		occi::Date computationdate(env, sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);

		int statNum = 0;
		IOUtils::convertString(statNum, stationNumber);

		stmt->setAutoCommit(true);

		unsigned int param = 1;

		stmt->setDate(param++, hazarddate);
		stmt->setString(param++, stationName);
		stmt->setNumber(param++, statNum);

		if (Hdata_ind[i].dewpt_def != -1)  stmt->setNumber(param++, Hdata[i].dewpt_def);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hoar_ind6 != -1)	stmt->setNumber(param++, Hdata[i].hoar_ind6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hoar_ind24 != -1) stmt->setNumber(param++, Hdata[i].hoar_ind24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wind_trans != -1) stmt->setNumber(param++, Hdata[i].wind_trans);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].hn3 != -1)       stmt->setNumber(param++, Hdata[i].hn3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hn6 != -1)       stmt->setNumber(param++, Hdata[i].hn6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hn12 != -1)      stmt->setNumber(param++, Hdata[i].hn12);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hn24 != -1)      stmt->setNumber(param++, Hdata[i].hn24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hn72 != -1)      stmt->setNumber(param++, Hdata[i].hn72);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hn72_24 != -1)   stmt->setNumber(param++, Hdata[i].hn72_24);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].hnw3 != -1)        stmt->setNumber(param++, Hdata[i].hnw3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hnw6 != -1)        stmt->setNumber(param++, Hdata[i].hnw6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hnw12 != -1)       stmt->setNumber(param++, Hdata[i].hnw12);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hnw24 != -1)       stmt->setNumber(param++, Hdata[i].hnw24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hnw72 != -1)       stmt->setNumber(param++, Hdata[i].hnw72);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].hoar_size != -1)  stmt->setNumber(param++, Hdata[i].hoar_size);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wind_trans24 != -1) stmt->setNumber(param++, Hdata[i].wind_trans24);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].stab_class1 != -1)  stmt->setNumber(param++, Hdata[i].stab_class1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_class2 != -1)  stmt->setNumber(param++, Hdata[i].stab_class2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index4);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height4);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index5);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height5);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].ch != -1)     stmt->setNumber(param++, Hdata[i].ch);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].crust != -1)  stmt->setNumber(param++, Hdata[i].crust);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].en_bal != -1) stmt->setNumber(param++, Hdata[i].en_bal);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].sw_net != -1)  stmt->setNumber(param++, Hdata[i].sw_net);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].t_top1 != -1)  stmt->setNumber(param++, Hdata[i].t_top1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].t_top2 != -1)  stmt->setNumber(param++, Hdata[i].t_top2);
		else stmt->setNull(param++, occi::OCCINUMBER);

		stmt->setNumber(param++, sn_version);
		stmt->setDate(param++, computationdate);

		if (Hdata_ind[i].swe != -1)      stmt->setNumber(param++, Hdata[i].swe);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].tot_lwc != -1)  stmt->setNumber(param++, Hdata[i].tot_lwc);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].runoff != -1)   stmt->setNumber(param++, Hdata[i].runoff);
		else stmt->setNull(param++, occi::OCCINUMBER);

		rows_inserted += stmt->executeUpdate(); // execute the statement stmt
		conn->terminateStatement(stmt);
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Inserted %d rows into %s!", rows_inserted, oracleDB.c_str());
}
