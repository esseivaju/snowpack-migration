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
/**
 * @file DataClasses.h
 * @version 11.08
 * This header file contains all the data structures needed for the 1d snowpack model
 */

#ifndef __DATACLASSES_H__
#define __DATACLASSES_H__

#include <snowpack/SnowpackConfig.h>

#include <snowpack/Constants.h>
#include <meteoio/MeteoIO.h>
#include <vector>
#include <sstream>

/// @brief The 3 different phases in the matrix
enum {
	SOLID,  ///< Solid
	LIQUID, ///< Liquid
	GAS,   ///< Gas
	NUMBER_OF_PHASES
};

/// @brief Optical Properties of snow
struct SNOW_OPTIC {
	double ggg;     ///< Asymmetry Parameter
	double exteff;  ///< Extinction Efficiency
	double ssa;     ///< Single Scattering Albedo
};

/// @brief Spectrum of incoming solar radiation
struct WL_STRUCT {
	double nm;      ///< Wavelength
	double perc;    ///< Percentage of Energy
};

/**
 * @brief ZwischenData contains "memory" information mainly for operational use\n
 * It is used to prepare some parameters of qr_Hdata. This data is read from and written to *.sno
 * or .haz files respectively
 */
class ZwischenData {
	public:
		ZwischenData(): hoar24(48, 0.0), drift24(48, 0.0), hn3(144, 0.0), hn24(144, 0.0) {};
		void reset();                ///< Sets all the values in the vectors to 0.0

		friend std::iostream& operator<<(std::iostream& os, const ZwischenData& data);
		friend std::iostream& operator>>(std::iostream& is, ZwischenData& data);

		std::vector<double> hoar24;  ///< Twenty-four hour hoar index every half-hour over one day 48
		std::vector<double> drift24; ///< Twenty-four hour hoar index every half-hour over one day 48
		std::vector<double> hn3;    ///< Three hour new snow heights every half-hour over three days 144
		std::vector<double> hn24;   ///< Twenty-four hour snow heights every half-hour over three days 144
};

/**
 * @brief CurrentMeteo is the class of interpolated meteo data for the current calculation time step \n
 * It contains some additional and very important derived parameters such as the roughness length or running mean values.
 */
class CurrentMeteo {
	public:
		CurrentMeteo();
		CurrentMeteo(const SnowpackConfig& i_cfg);
		void reset(const SnowpackConfig& i_cfg);
		void setMeasTempParameters(const mio::MeteoData& md);
		size_t getNumberMeasTemperatures() const;
		size_t getNumberFixedRates() const;
		size_t getMaxNumberMeasTemperatures() const;
		void getFixedPositions(std::vector<double>& positions) const;
		size_t getNumberFixedPositions() const;
		void copySnowTemperatures(const mio::MeteoData& md, const unsigned int& current_slope);
		void copySolutes(const mio::MeteoData& md, const size_t& i_number_of_solutes);

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const CurrentMeteo& data);
		friend std::iostream& operator>>(std::iostream& is, CurrentMeteo& data);

		mio::Date date;  ///< Date of current meteo data
		double ta;       ///< Air temperature (K)
		double rh;       ///< Relative humidity (% or 1)
		double rh_avg;   ///< Running mean of relative humidity (1)
		double vw;       ///< Wind velocity at snow station (m s-1)
		double vw_avg;   ///< Running mean of wind velocity at snow station (m s-1)
		double vw_max;   ///< Maximum wind velocity at snow station (m s-1)
		double dw;       ///< Wind direction at snow station (deg)
		double vw_drift; ///< Wind velocity for blowing and drifting snow (operational: wind ridge station)
		double dw_drift; ///< Wind direction of blowing and drifting snow (operational: wind ridge station)
		double ustar;    ///< The friction velocity (m s-1) computed in mt_MicroMet() and also used later for the MeteoHeat fluxes
		double z0;       ///< The roughness length computed in SnowDrift and also used later for the MeteoHeat fluxes (m)
		double psi_s;    ///< Stability correction for scalar heat fluxes
		double iswr;     ///< Incoming SHORTWAVE radiation (W m-2)
		double rswr;     ///< Reflected SHORTWAVE radiation (W m-2) divide this value by the ALBEDO to get iswr
		double mAlbedo;  ///< Measured snow albedo
		double diff;     ///< Diffuse radiation from the sky (W m-2)
		double dir_h;    ///< Horizontal direct radiation from the sky (W m-2)
		double elev;     ///< Solar elevation to be used in Canopy.c (rad) => see also
		double ea;       ///< Atmospheric emissivity (1)
		double tss;      ///< Snow surface temperature (K)
		double tss_a12h; ///< Snow surface temperature averaged over past 12 hours (K)
		double tss_a24h; ///< Snow surface temperature averaged over past 24 hours (K)
		double ts0;      ///< Bottom temperatures of snow/soil pack (K)
		double hnw;      ///< The water equivalent of snowfall in mm w.e. (kg m-2) per CALCULATION_STEP_LENGTH
		double hs;       ///< The measured height of snow (m)
		double hs_a3h;   ///< Snow depth averaged over 3 past hours
		double hs_rate;  ///< The rate of change in snow depth (m h-1)
		double adv_heat; ///< Advective heat to inject in the soil (if ADVECTIVE_HEAT and related parameters set to true)

		std::vector<double> ts;    ///< Measured snow or/and soil temperatures (K)
		std::vector<double> zv_ts; ///< Positions of all measured snow or/and soil temperatures (m)
		std::vector<double> conc;  ///< Solute concentrations in precipitation
		double rho_hn;             ///< Measured new snow density (kg m-3)

	private:
		size_t getNumberMeasTemperatures(const mio::MeteoData& md);

		std::vector<double> fixedPositions; ///< Positions of fixed snow/soil temperatures (m)
		double minDepthSubsurf;             ///< Sensor must be covered by minDepthSubsurf (m) to be output
		size_t maxNumberMeasTemperatures;   ///< Max allowed number of measured snow/soil temperatures, depending on variant
		size_t numberMeasTemperatures;      ///< Number of measured snow/soil temperatures
		size_t numberFixedRates;
};

/// @brief The 3 mathematical fields that can be solved
enum SN_FIELD{
	TEMPERATURE, ///< Temperature (phase change and metamorphism)
	SEEPAGE,     ///< Water transport
	SETTLEMENT,  ///< Creep displacements
	N_SN_FIELDS
};

/// @brief The 4 different components of the porous matrix
enum {
	SOIL,  ///< Soil
	ICE,   ///< Ice
	WATER, ///< Water
	AIR,   ///< Air
	N_COMPONENTS
};

/// @brief The different soil fields (PERMAFROST)
enum SN_SOIL_DATA{
	SOIL_RHO, ///< Bulk density of dry soil component (without accounting for porosity)
	SOIL_K,   ///< Bulk conductivity of dry soil component
	SOIL_C,   ///< Specific heat of dry soil component
	N_SOIL_FIELDS
};

/**
 * @brief Parameters of the different layers of the snowpack \n
 * It is used as a pointer (array) within the SSdata (profile) data structure.
 */
class LayerData {
	public:
		LayerData();

		friend std::iostream& operator<<(std::iostream& os, const LayerData& data);
		friend std::iostream& operator>>(std::iostream& is, LayerData& data);

		mio::Date depositionDate;   ///< Date of deposition (mainly used for snow layers)
		double hl;                  ///< The height of the layer in m
		size_t ne;                  ///< Number of finite elements in the the layer (hl/ne defines elm. size)
		double tl;                  ///< Temperature at the top of the layer in K or degC
		double phiSoil;             ///< Volumetric soil content in %
		double phiIce;              ///< Volumetric ice content in %
		double phiWater;            ///< Volumetric water content in %
		double phiVoids;            ///< Volumetric void content in %
		std::vector<double> cSoil;  ///< Solute concentrations in Soil
		std::vector<double> cIce;   ///< Solute concentrations in Ice
		std::vector<double> cWater; ///< Solute concentrations in Water
		std::vector<double> cVoids; ///< Solute concentrations in Air
		double SoilRho;             ///< Density of soil in kg m-3
		double SoilK;               ///< Conductivity of soil
		double SoilC;               ///< Heat Capacity of soil
		double rg;                  ///< Micro-structure : Grainsize in mm
		double sp;                  ///< Micro-structure : Sphericity
		double dd;                  ///< Micro-structure : Dendricity
		double rb;                  ///< Micro-structure : Bond Radius in mm
		unsigned short int mk;      ///< Micro-structure : Marker
		double hr;                  ///< Surface hoar Mass in kg m-2
		double CDot;                ///< Stress rate (Pa s-1), that is the LAST overload change rate
		double metamo;              ///< keep track of metamorphism
};

/**
 * @brief SN_SNOWSOIL_DATA includes all important station parameters as well as LayerData \n
 * This data structure will have to be replaced by something a little more complicated soon ???
 * For now it is simply an efficient way of creating a snowpack to investigate.
 */
class SN_SNOWSOIL_DATA {
	public:
		SN_SNOWSOIL_DATA() : meta(), profileDate(), nN(0), Height(0.),
                     nLayers(0), Ldata(), HS_last(0.), Albedo(0.), SoilAlb(0.), BareSoil_z0(0.),
                     Canopy_Height(0.), Canopy_LAI(0.), Canopy_Direct_Throughfall(0.),
                     WindScalingFactor(1.), ErosionLevel(0), TimeCountDeltaHS(0.)
		{
			Ldata.clear();
		}

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const SN_SNOWSOIL_DATA& data);
		friend std::iostream& operator>>(std::iostream& is, SN_SNOWSOIL_DATA& data);

		mio::StationData meta;            ///< Station meta data
		mio::Date profileDate;            ///< Date of profile
		size_t nN;                        ///< Total number of FE nodes
		double Height;                    ///< Total height of snowpack in m (sum of the layer heights)
		size_t nLayers;                   ///< Total number of snowpack layers
		std::vector<LayerData> Ldata;     ///< contains all the information required to construct the Xdata
		double HS_last;                   ///< Last checked measured Snow Height
		double Albedo;                    ///< Snow albedo
		double SoilAlb;                   ///< Soil albedo; default 0.2
		double BareSoil_z0;               ///< Bare soil roughness in m; default 0.02 m
		double Canopy_Height;             ///< Canopy Height in m
		double Canopy_LAI;                ///< Canopy Leaf Area Index in m2 m-2
		double Canopy_Direct_Throughfall; ///< Direct throughfall [fraction of precipitation]
		double WindScalingFactor;         ///< Local scaling factor for wind at drift station
		int    ErosionLevel;              ///< Erosion Level in operational mode (flat field virtual erosion)
		double TimeCountDeltaHS;          ///< Time counter tracking erroneous settlement in operational mode
};

/**
 * @brief ELEMENT DATA used as a pointer in the SnowStation structure
 * NOTE on M below: this is the mass of an element that is neither changed by phase changes nor densification. \n
 * It is set in the data initialization and used to compute the stress field.
 * It can ONLY be changed by the WATER TRANSPORT or SURFACE SUBLIMATION or WIND TRANSPORT routines.
 */
class ElementData {
	public:
		ElementData();

		bool checkVolContent() const;
		void heatCapacity();
		double coldContent() const;
		double extinction() const;
		void opticalEquivalentGrainSize();
		void snowResidualWaterContent();
		static double snowResidualWaterContent(const double& theta_i);
		double soilFieldCapacity() const;

		double snowElasticity() const;
		double neckStressEnhancement() const;
		double concaveNeckRadius() const;
		double neckLength() const;
		double neck2VolumetricStrain() const;

		void snowType();
		static unsigned short int snowType(const double& dendricity, const double& sphericity, const double& grain_dia, const size_t& marker,
                        const double& theta_w, const double& res_wat_cont);

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const ElementData& data);
		friend std::iostream& operator>>(std::iostream& is, ElementData& data);

		mio::Date depositionDate;  ///< Date of deposition
		double L0, L;              ///< Original and present element length (m)
		double Te;                 ///< mean element temperature (K)
		double gradT;              ///< temperature gradient over element (K m-1)
		double melting_tk;	   ///< melt temperature of layer (principally initialized as 0 degC, but enables possibility for freezing point depression)
		double freezing_tk;	   ///< freezing temperature of layer (principally initialized as 0 degC, but enables possibility for freezing point depression)
		std::vector<double> theta; ///< volumetric contents: SOIL, ICE, WATER, AIR (1)
		mio::Array2D<double> conc; ///< Concentration for chemical constituents in (kg m-3)
		std::vector<double> k;     ///< For example, heat conductivity of TEMPERATURE field (W m-1 K-1)
		//   Stored in order to visualize constitutive laws
		//   Will be used for creep field hydraulic conductivity in m3 s kg-1
		std::vector<double> c;     ///< For example, specific heat of TEMPERATURE field (J kg-1 K-1)
		//   Will also be used for creep specific snow water capacity  in m3 J-1
		std::vector<double> soil;  ///< Contains the heat conductivity, capacity and dry density of the soil (solid, non-ice)  component phase
		double Rho;                ///< mean element density (or BULK density; kg m-3), that is, rho=M/V=sum( theta(i)*rho(i) )
		double M;                  ///< the total mass of the element (kg m-2)
		double sw_abs;             ///< total absorbed shortwave radiation by the element (W m-2)
		// Snow Metamorphism Data
		double rg;                 ///< grain radius (mm)
		double dd;                 ///< snow dendricity: 0 = none, 1 = newsnow
		double sp;                 ///< sphericity: 1 = round, 0 = angular
		double ogs;                ///< optical equivalent grain size (mm)
		double rb;                 ///< grain bond radius (mm)
		double N3;                 ///< grain Coordination number (1)
		size_t mk;                 ///< grain marker (history dependent)
		unsigned short int type;   ///< grain class
		double metamo;             ///< keep track of metamorphism
		double dth_w;              ///< Subsurface Melting & Freezing Data: change of water content
		double res_wat_cont;       ///< Residual water content
		double Qmf;                ///< Subsurface Melting & Freezing Data: change of energy due to phase changes (melt-freeze)
		double dE, E, Ee, Ev;      ///< Total element strain (GREEN'S strains -- TOTAL LAGRANGIAN FORMULATION.
		double EDot, EvDot;        ///< Total Strain Rate, elastic and viscous, respectively (s-1) (Simply, E/sn_dt)
		double S;                  ///< Total Element Stress (Pa), S being the energy conjugate stress
		double C;                  ///< Total Element Stress (Pa), C being the real or the Cauchy stress, which is output
		double CDot;               ///< Stress rate (Pa s-1), that is the overload change rate
		double ps2rb;              ///< proportion of grain bond growth due to pressure sintering (1)
		double s_strength;         ///< Parameterized snow shear strength (kPa)
		double hard;               ///< Parameterized hand hardness (1)
		double S_dr;               ///< Stability Index based on deformation rate (Direct Action Avalanching)
		double theta_r;            ///< Residual water content of previous time step (m^3/m^3), used exclusively for solving Richards equation in snow
		//NIED (H. Hirashima)
		double dhf;
};

/// @brief NODAL DATA used as a pointer in the SnowStation structure
class NodeData {
	public:
		NodeData() : z(0.), u(0.), f(0.), udot(0.), T(0.), S_n(0.), S_s(0.), ssi(6.), hoar(0.),
		             dhf(0.), S_dhf(0.), Sigdhf(0.) {} //HACK: set ssi to max_stability!

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const NodeData& data);
		friend std::iostream& operator>>(std::iostream& is, NodeData& data);

		double z;    ///< nodal height from ground in m
		double u;    ///< creep displacements in m
		double f;    ///< reaction or unbalanced forces (CREEP)
		double udot; ///< downward creep velocity in m s-1
		double T;    ///< nodal temperature in K
		double S_n;  ///< Stability Index for natural avalanches
		double S_s;  ///< Stability Index for skier triggered avalanches
		double ssi;  ///< Structural Stability Index
		double hoar; ///< Mass of surface hoar collected while node was exposed to surface

		//NIED (H. Hirashima)
		double dhf;
		double S_dhf;
		double Sigdhf;
};

/**
 * @brief Canopy data used as a pointer in the SnowStation structure
 * -# INSTANTANEOUS VARIABLES
 * 	-# Canopy "state" variables, and some auxiliaries
 * 	-# Properties which could be given here or as a parameter field
 * 	-# Aerodynamic resistances
 * -# CUMULATED/AVERAGE VARIABLES \n
 *    cumulated between time series output timesteps - these variables can be moved or mirrored in a canopy "surface" data structure
 * 	-# Albedo and similar auxiliaries
 * 	-# Radiation fluxes
 * 	-# Canopy turbulent heat fluxes
 * 	-# Canopy evaporative fluxes
 * 	-# Canopy mass fluxes
 */
class CanopyData {
	public:
		CanopyData() : storage(0.), temp(0.), sigf(0.), ec(0.), lai(0.), z0m(0.), z0h(0.), zdispl(0.),
		     height(0.), direct_throughfall(0.), ra(0.), rc(0.), rs(0.), rstransp(0.), canopyalb(0.),
		     totalalb(0.), wetfraction(0.), intcapacity(0.), rswrac(0.), iswrac(0.), rswrbc(0.),
		     iswrbc(0.), ilwrac(0.), rlwrac(0.), ilwrbc(0.), rlwrbc(0.), rsnet(0.), rlnet(0.),
		     sensible(0.), latent(0.), latentcorr(0.), transp(0.), intevap(0.),
		     interception(0.), throughfall(0.), snowunload(0.) {}

		void reset(const bool& cumsum_mass);
		void initializeSurfaceExchangeData();

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const CanopyData& data);
		friend std::iostream& operator>>(std::iostream& is, CanopyData& data);

		// Aa
		double storage;  ///< intercepted water (mm or kg m-2)
		double temp;	   ///< temperature (K)
		double sigf;	   ///< radiation transmissivity (1)
		double ec;       ///< longwave emissivity (1)
		// Ab
		double lai;
		double z0m;
		double z0h;
		double zdispl;
		double height;
		double direct_throughfall;
		// Ac
		double ra;          ///< from canopy air to reference height
		double rc;          ///< from canopy to canopy air
		double rs;          ///< from subsurface to canpopy air
		double rstransp;    ///< stomatal surface resistance for transpiration
		// Ba
		double canopyalb;   ///< canopy albedo [-]
		double totalalb;    ///< total albedo above canopy and snow/soil surface [-]
		double wetfraction; ///< fraction of canopy covered by interception [-]
		double intcapacity; ///< maximum interception storage [mm]
		// Bb
		double rswrac;      ///< upward shortwave above canopy
		double iswrac;	    ///< downward shortwave radiation above canopy
		double rswrbc;      ///< upward shortwave below canopy
		double iswrbc;      ///< downward shortwave radiation below canopy
		double ilwrac;      ///< downward longwave radiation ABOVE canopy
		double rlwrac;      ///< upward longwave radiation ABOVE canopy
		double ilwrbc;      ///< downward longwave radiation BELOW canopy
		double rlwrbc;      ///< upward longwave radiation BELOW canopy
		double rsnet;       ///< net shortwave radiation
		double rlnet;       ///< net longwave radiation
		// Bc
		double sensible;
		double latent;
		double latentcorr;
		// Bd
		double transp;
		double intevap;
		// Be
		double interception;
		double throughfall;
		double snowunload;
};

/**
 * @brief Station data including all information on snowpack layers (elements and nodes) and on canopy \n
 * This is the PRIMARY data structure of the SNOWPACK program \n
 * It is used extensively not only during the finite element solution but also to control
 * the post-processing writes. It is initialized from SN_SNOWSOIL_DATA (at present).
 */
class SnowStation {
	public:
		explicit SnowStation(const bool& i_useCanopyModel=true, const bool& i_useSoilLayers=true);
		SnowStation(const SnowStation& c);

		~SnowStation();
		SnowStation& operator=(const SnowStation&); ///<Assignement operator

		void initialize(const SN_SNOWSOIL_DATA& SSdata, const size_t& i_sector);
		void resize(const size_t& number_of_elements);

		void reduceNumberOfElements(const size_t& rnE);
		void combineElements(const size_t& number_top_elements);
		static bool combineCondition(const ElementData& Edata0, const ElementData& Edata1);
		static void mergeElements(ElementData& Edata0, const ElementData& Edata1, const bool& merge, const bool& topElement);

		void compSnowpackMasses();
		void compSnowpackInternalEnergyChange(const double& sn_dt);
		void compSoilInternalEnergyChange(const double& sn_dt);
		double getLiquidWaterIndex() const;
		double getModelledTemperature(const double& z) const;

		size_t getNumberOfElements() const;
		size_t getNumberOfNodes() const;
		bool isGlacier(const bool& hydro=false) const;
		bool hasSoilLayers() const;

		size_t find_tag(const size_t& tag) const;

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const SnowStation& data);
		friend std::iostream& operator>>(std::iostream& is, SnowStation& data);

		mio::StationData meta;      ///< Station meta data
		double cos_sl;              ///< Cosinus of slope angle, initialized once!
		size_t sector;              ///< current slope sector of width 360./MAX(1, nSlopes-1)

		CanopyData Cdata;           ///< Pointer to canopy data
		double pAlbedo;             ///< Parameterized snow albedo
		double Albedo;              ///< Snow albedo used by the model
		double SoilAlb;             ///< Soil albedo
		double BareSoil_z0;         ///< Bare soil roughness in m
		size_t SoilNode;            ///< The top soil node, 0 in case of SNP_SOIL == 0
		double Ground;              ///< The ground height -- meaning the height of the top soil node
		double cH;                  ///< The CALCULATED snowpack height, including soil depth if SNP_SOIL == 1
		double mH;                  ///< The ENFORCED snowpack height, including soil depth if SNP_SOIL == 1
		double mass_sum;            ///< Total mass summing mass of elements
		double swe;                 ///< Total mass summing snow water equivalent of elements
		double lwc_sum;             ///< Total liquid water in snowpack
		double hn;                  ///< Depth of new snow to be used on slopes
		double rho_hn;              ///< Density of new snow to be used on slopes
		size_t ErosionLevel;        ///< Element where snow erosion stopped previously for the drift index
		double ErosionMass;         ///< Eroded mass either real or virtually (storage if less than one element)
		int S_class1;               ///< Stability class based on hand hardness, grain class ...
		int S_class2;               ///< Stability class based on hand hardness, grain class ...
		double S_d;                 ///< Minimum Direct Action Stability Index  ...
		double z_S_d;               ///< Depth of Minimum Direct Action Stability
		double S_n;                 ///< Minimum Natural Stability Index
		double z_S_n;               ///< Depth of Minimum Natural Stability
		double S_s;                 ///< Minimum Alternative Stability Index (ASI; skier stability)
		double z_S_s;               ///< Depth of Minimum ASI
		double S_4;                 ///< stab_index4
		double z_S_4;               ///< Depth of stab_index4
		double S_5;                 ///< stab_index5
		double z_S_5;               ///< Depth of stab_index5
		std::vector<NodeData> Ndata;    ///< pointer to nodal data array (e.g. T, z, u, etc..)
		std::vector<ElementData> Edata; ///< pointer to element data array (e.g. Te, L, Rho, etc..)
		void *Kt;                   ///< Pointer to pseudo-conductivity and stiffnes matrix
		size_t tag_low;             ///< Lowest tag to dump, 0 means no tags at all
		double ColdContent;         ///< Cold content of snowpack (J m-2)
		double ColdContentSoil;     ///< Cold content of soil (J m-2)
		double dIntEnergy;          ///< Internal energy change of snowpack (J m-2)
		double dIntEnergySoil;      ///< Internal energy change of soil (J m-2)
		double meltFreezeEnergy;    ///< Melt freeze part of internal energy change of snowpack (J m-2)
		double meltFreezeEnergySoil;///< Melt freeze part of internal energy change of soil (J m-2)
		double ReSolver_dt;         ///< Last used RE time step in the previous SNOWPACK time step
		char SubSurfaceMelt;        ///< Subsurface melting flag ( yes/no ) for exposition
		char SubSurfaceFrze;        ///< Subsurface refreezing flag ( yes/no ) for exposition
		bool windward;              ///< True for windward (luv) slope

		static const double comb_thresh_l, comb_thresh_ice, comb_thresh_water;
		static const double comb_thresh_dd, comb_thresh_sp, comb_thresh_rg;
		static const double thresh_moist_snow, thresh_moist_soil;
		static const size_t number_top_elements;
		static unsigned short number_of_solutes;  ///< The model treats that number of solutes

	private:
		size_t nNodes;                      ///< Actual number of nodes; different for each exposition
		size_t nElems;                      ///< Actual number of elements (nElems=nNodes-1)
		bool useCanopyModel, useSoilLayers; ///< The model includes soil layers
};

/**
* @brief BoundCond is used to set Neumann boundary conditions
*/
class BoundCond {

	public:
		BoundCond() : lw_out(0.), lw_net(0.), qs(0.), ql(0.), qr(0.), qg(Constants::undefined) {};

		double lw_out;  ///< outgoing longwave radiation
		double lw_net;  ///< net longwave radiation
		double qs;      ///< sensible heat
		double ql;      ///< latent heat
		double qr;      ///< rain energy
		double qg;      ///< geothermal heat flux or heat flux at lower boundary
};

/**
* @name Surface data
* @note Some of the most important results of the simulation are contained in these data structures
*/
//@{
class SurfaceFluxes {
	public:
		/**
		 * @brief The different types of mass fluxes:
		 * Mass fluxes in kg m-2 \n
		 * Rates in kg m-2 h-1 (MS_HNW, MS_RAIN and MS_WIND)
		 */
		enum SN_MASS_CHANGES {
		MS_TOTALMASS,      ///< This of course is the total mass of the snowpack at the present time
		MS_SWE,            ///< This too, of course, but summing rho*L
		MS_WATER,          ///< The total amount of water in the snowpack at the present time
		MS_HNW,            ///< Solid precipitation rate
		MS_RAIN,           ///< Rain rate
		MS_WIND,           ///< Mass loss rate due to wind erosion
		MS_EVAPORATION,    ///< The mass loss or gain of the top element due to water evaporating
		MS_SUBLIMATION,    ///< The mass loss or gain of the top element due to snow (ice) sublimating
		MS_SNOWPACK_RUNOFF,///< The total mass loss of snowpack due to water transport (virtual lysimeter)
		MS_SOIL_RUNOFF,    ///< Equivalent to MS_SNOWPACK_RUNOFF but at bottom soil node
		N_MASS_CHANGES     ///< Total number of different mass change types
		};

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const SurfaceFluxes& data);
		friend std::iostream& operator>>(std::iostream& is, SurfaceFluxes& data);

		SurfaceFluxes();

		void reset(const bool& cumsum_mass);
		void compSnowSoilHeatFlux(const SnowStation& Xdata);
		void collectSurfaceFluxes(const BoundCond& Bdata, SnowStation& Xdata, const CurrentMeteo& Mdata);

		/**
		 * @brief Energy fluxes:
		 * Energy change of snowpack in kJ m-2 (dIntEnergy)\n
		 * Fluxes in W m-2
		 */
		double lw_in;      ///< incoming longwave radiation
		double lw_out;     ///< outgoing longwave radiation
		double lw_net;     ///< net longwave radiation
		double qs;         ///< sensible heat
		double ql;         ///< latent heat
		double hoar;       ///< mass of surface hoar formed or sublimated
		double qr;         ///< rain energy
		double qg;         ///< geothermal heat flux or heat flux at lower boundary
		double qg0;        ///< ground heat flux at soil-snow interface
		double sw_hor;     ///< incoming global shortwave radiation on horizontal surface
		double sw_in;      ///< incoming global shortwave radiation; on slopes: projected
		double sw_out;     ///< reflected shortwave radiation
		double qw;         ///< net shortwave radiation at the surface (absorbed within the snowpack)
		double sw_dir;     ///< incoming direct shortwave radiation; on slopes: projected
		double sw_diff;    ///< incoming diffuse shortwave radiation
		double pAlbedo;    ///< parameterized Albedo (used for OUTPUT only)
		double mAlbedo;    ///< measured Albedo (used for OUTPUT only)
		double dIntEnergy;           ///< Internal energy change in J m-2 in snowpack (used for OUTPUT only)
		double dIntEnergySoil;       ///< Internal energy change in J m-2 in soil (used for OUTPUT only)
		double meltFreezeEnergy;     ///< Melt freeze part of internal energy change in J m-2 in snowpack (used for OUTPUT only)
		double meltFreezeEnergySoil; ///< Melt freeze part of internal energy change in J m-2 in soil (used for OUTPUT only)

		/// @brief Other surface data:
		double drift;      ///< the surface flux of drifting snow in kg m-1 s-1
		std::vector<double> mass; ///< Total mass of snowpack PLUS different amounts of total mass change, sublimation, runoff, erosion, etc. Basically the mass which crosses the surface
		std::vector<double> load; ///< Total load (kg m-2) in water runoff from solutes like nitrate
		double dhs_corr;   ///< operational mode only: snow depth correction in case of squeezing or blow-up (m)
		double cRho_hn;    ///< Computed new snow density (kg m-3)
		double mRho_hn;    ///< Measured new snow density (kg m-3)

};
//@}

/// @brief Defines structure for snow profile layers
class SnowProfileLayer {
	public:
		SnowProfileLayer();

		void average(const double& w1, const double& w2, const SnowProfileLayer& Pdata);
		static std::vector<SnowProfileLayer> generateProfile(const mio::Date& dateOfProfile, const SnowStation& Xdata, const double hoar_density_surf, const double hoar_min_size_surf);

		// Profile meta data
		mio::Date profileDate; ///< Date of profile
		std::string stationname;
		size_t  loc_for_snow;
		size_t  loc_for_wind;

		mio::Date depositionDate;   ///< Date of deposition (mainly used for snow layers)
		double height;         ///< 0 to 1000      (cm)
		double rho;            ///< 0 to 1000      (kg m-3)
		double T;              ///< -50 to 50, snow temperature at top of layer (degC)
		double gradT;          ///< -1000 to 1000, temperature gradient across layer (K m-1)
		double v_strain_rate;  ///< 0 to 1.0e-5, viscous strain rate (s-1)
		double theta_i;        ///< 0 to 1, volume fraction of ice (-)
		double theta_w;        ///< 0 to 1, volume fraction of water (-)
		double theta_a;        ///< 0 to 1, volume fraction of air (-)
		double grain_size;     ///< 0 to 100       (mm)
		double bond_size;      ///< 0 to 100       (mm)
		double dendricity;     ///< 0 to 1         (-)
		double sphericity;     ///< 0 to 1         (-)
		double ogs;            ///< 0 to 100, optical equivalent grain size (mm)
		double coordin_num;    ///< 0 to 10        (-)
		size_t marker;         ///< 0 to 999       (-)
		short unsigned int type; ///< 0 to 999     (-)
		double hard;           ///< 0. to 5.       (-)

	private:
		void generateLayer(const ElementData& Edata, const NodeData& Ndata);
		void generateLayer(const ElementData& Edata, const NodeData& Ndata,
                           const mio::Date& dateOfProfile, const double hoar_density_surf);
};

/// @brief class to collect the information about the current simulation (version, date)
class RunInfo {
	public:
		RunInfo();

		const std::string version;   ///< SNOWPACK version
		const mio::Date computation_date; ///< Date of computation
		const std::string compilation_date; ///< Date of compilation
		const std::string user; ///< logname of the user running the simulation

	private:
		static mio::Date getRunDate();
		static std::string getCompilationDate();
};

/// Structure of double values for output to SDB
struct ProcessDat {
	ProcessDat() : date(), nHz(0), stat_abbrev(), loc_for_snow(0), loc_for_wind(0),
	               ch(0.), swe(0.), tot_lwc(0.), runoff(0.), dewpt_def(0.), hoar_size(0.), hoar_ind6(0.), hoar_ind24(0.),
	               wind_trans(0.), wind_trans24(0.),
	               hn_half_hour(0.), hn3(0.), hn6(0.), hn12(0.), hn24(0.), hn72(0.), hn72_24(0.),
	               hnw_half_hour(0.), hnw3(0.), hnw6(0.), hnw12(0.), hnw24(0.), hnw72(0.),
	               stab_class1(0), stab_class2(0),
	               stab_index1(0.), stab_height1(0.), stab_index2(0.), stab_height2(0.), stab_index3(0.), stab_height3(0.), stab_index4(0.),stab_height4(0.), stab_index5(0.), stab_height5(0.),
	               crust(0.), en_bal(0.), sw_net(0.), t_top1(0.), t_top2(0.), lwi_N(0.), lwi_S(0.),
	               dhs_corr(0.), mass_corr(0.)
	{}

	mio::Date date;        ///< Process date
	int nHz;               ///< Number of hazard steps
	char stat_abbrev[16];
	int  loc_for_snow;
	int  loc_for_wind;
	// Data
	double ch;             ///< height of snow HS (cm)
	double swe;            ///< snow water equivalent SWE (kg m-2)
	double tot_lwc;        ///< total liquid water content (kg m-2)
	double runoff;         ///< runoff (kg m-2)
	double dewpt_def;      ///< dew point deficit (degC)
	double hoar_size;      ///< 24 h surface hoar size (mm)
	double hoar_ind6;      ///< 6 h surface hoar index (kg m-2)
	double hoar_ind24;     ///< 24 h surface hoar index (kg m-2)
	double wind_trans;     ///< 6 h drifting snow index (cm)
	double wind_trans24;   ///< 24 h drifting snow index (cm)
	double hn_half_hour;   ///< half_hour depth of snowfall (cm)
	double hn3;            ///< 3 h depth of snowfall (cm)
	double hn6;            ///< 6 h depth of snowfall (cm)
	double hn12;           ///< 12 h depth of snowfall (cm)
	double hn24;           ///< 24 depth of snowfall (cm)
	double hn72;           ///< 72 depth of snowfall (cm)
	double hn72_24;        ///< 3 d sum of 24 h depth of snowfall (cm)
	double hnw_half_hour;  ///< half_hour new snow water equivalent (kg m-2)
	double hnw3;           ///< 3 h new snow water equivalent (kg m-2)
	double hnw6;           ///< 6 h new snow water equivalent (kg m-2)
	double hnw12;          ///< 12 h new snow water equivalent (kg m-2)
	double hnw24;          ///< 24 h new snow water equivalent (kg m-2)
	double hnw72;          ///< 72 h new snow water equivalent (kg m-2)
	int stab_class1;       ///< stability classes 1,3,5
	int stab_class2;       ///< profile type 0..10
	double stab_index1;    ///< deformation index Sdef
	double stab_height1;   ///< depth of stab_index1 (cm)
	double stab_index2;    ///< natural stability index Sn38
	double stab_height2;   ///< depth of stab_index2 (cm)
	double stab_index3;    ///< skier stability index Sk38
	double stab_height3;   ///< depth of stab_index3 (cm)
	double stab_index4;    ///< structural stability index SSI
	double stab_height4;   ///< depth of stab_index4 (cm)
	double stab_index5;    ///< none
	double stab_height5;   ///< depth of stab_index5 (cm)
	// Special parameters
	double crust;          ///< height of melt-freeze crust on southern slope (cm)
	double en_bal;         ///< internal energy change (kJ m-2)
	double sw_net;         ///< surface energy input (kJ m-2)
	double t_top1, t_top2; ///< snow temperatures at depth 1 & 2, respectively (degC)
	double lwi_N, lwi_S;   ///< liquid water index for northerly and southerly slopes, respectively.
	// Control parameters
	double dhs_corr;  ///< snow depth correction in case of squezzing or blow-up (cm)
	double mass_corr; ///< mass correction from either forced erosion and squeezing (neg) or blowing up (pos) (cm)
};

struct ProcessInd {
	ProcessInd() : stat_abbrev(0), loc_for_snow(0), loc_for_wind(0),
	               ch(0), swe(0), tot_lwc(0), runoff(0), dewpt_def(0),
	               hoar_size(0), hoar_ind6(0), hoar_ind24(0),
	               wind_trans(0), wind_trans24(0),
	               hn3(0), hn6(0), hn12(0), hn24(0), hn72(0), hn72_24(0), hnw3(0), hnw6(0), hnw12(0), hnw24(0), hnw72(0),
	               stab_class1(0), stab_class2(0),
	               stab_index1(0), stab_height1(0), stab_index2(0), stab_height2(0), stab_index3(0), stab_height3(0), stab_index4(0), stab_height4(0), stab_index5(0), stab_height5(0),
	               crust(0), en_bal(0), sw_net(0), t_top1(0), t_top2(0), lwi_N(0), lwi_S(0)
	{}

	short stat_abbrev;
	short loc_for_snow;
	short loc_for_wind;
	// Data
	short ch;
	short swe;
	short tot_lwc;
	short runoff;
	short dewpt_def;
	short hoar_size;
	short hoar_ind6, hoar_ind24;
	short wind_trans, wind_trans24;
	short hn3, hn6, hn12, hn24, hn72;
	short hn72_24;
	short hnw3, hnw6, hnw12, hnw24, hnw72;
	short stab_class1, stab_class2;
	short stab_index1, stab_height1;
	short stab_index2, stab_height2;
	short stab_index3, stab_height3;
	short stab_index4, stab_height4;
	short stab_index5, stab_height5;
	short crust;
	short en_bal;
	short sw_net;
	short t_top1, t_top2;
	short lwi_N, lwi_S;
};

/// @brief Class for recording reference properties of tagged elements
class Tag {
	public:
		Tag();

		void compute_properties(const ElementData& Edata);
		void reposition_tag(const bool& useSoilLayers, const double& z, SnowStation& Xdata);

		static const bool metamo_expl; ///< set while using the explicit metamorphism model

		std::string label;             ///< Label for output file header
		mio::Date date;                ///< date at which to start tagging
		//char label[MAX_STRING_LENGTH]; ///< Label for output file header
		//double JulianDate;             ///< Julian date at which to start tagging

		// Repositioning
		size_t elem;                      ///< Index of tagged element
		double previous_depth;         ///< Last position of corresponding fixed rate sensor perpendicular to slope (m)
		// Viscosity
		double etaNS;                  ///< New snow viscosity according to M. Lehning
		double etaMSU;                 ///< Snow viscosity (Montana model)
		// Metamorphism
		double ML2L;                   ///< layer to layer flux
		double lp;                     ///< lattice constant
};

class TaggingData {
	public:
		TaggingData(const double& i_calculation_step_length);
		void resize(size_t i_size);
		void update_tags(const CurrentMeteo&  Mdata, SnowStation& Xdata);

		bool useSoilLayers, surface_write;
		double calculation_step_length;
		size_t tag_low, tag_top, repos_low, repos_top;
		std::vector<Tag> tags;

	private:
		size_t number_tags;
};

#endif
