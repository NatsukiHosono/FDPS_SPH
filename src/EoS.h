#pragma once

namespace EoS{
	//////////////////
	//abstract class
	//////////////////
	template <typename type> class EoS_t{
		public:
		EoS_t(){
			return ;
		}
		/* virtual */ ~EoS_t(){
			return ;
		}
		virtual type Pressure  (const type dens, const type eng) const = 0;
		virtual type SoundSpeed(const type dens, const type eng) const = 0;
	};
	//////////////////
	//EoSs
	//////////////////
	template <typename type> class IdealGas : public EoS_t<type>{
		const type hcr;//heat capacity ratio;
		public:
		IdealGas(const type _hcr) : hcr(_hcr){
		}
		inline type Pressure(const type dens, const type eng) const{
			return (hcr - 1.0) * dens * eng;
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			return sqrt(hcr * (hcr - 1.0) * eng);
		}
		inline type HeatCapacityRatio() const{
			return hcr;
		}
	};
<<<<<<< HEAD
=======

>>>>>>> a4dbabc8242229d27592959090ddc38d4abd65dd
    //Tillotson equation "is applicable to the prediction of the shock and release of materials undergoing hypervelocity impacts."
	template <typename type> class Tillotson : public EoS_t<type>{
		type rho0, a, b, A, B, u0, alpha, beta, uiv, ucv;
		inline type P_co(const type dens, const type eng) const{
			const type eta = dens / rho0;
			const type mu  = eta - 1.0;
			return (a + b / (eng / u0 / eta / eta + 1.0)) * dens * eng + A * mu + B * mu * mu;
		}
		inline type P_ex(const type dens, const type eng) const{
			const type eta = dens / rho0;
			const type mu  = eta - 1.0;
			return a * dens * eng + (b * dens * eng / (eng / u0 / eta / eta + 1.0) + A * mu * exp(- alpha * (1.0 / eta - 1.0))) * exp(- beta * (1.0 / eta - 1.0) * (1.0 / eta - 1.0));
		}
		inline type dPdrho(const type rho, const type u) const{
			const type drho = 0.0001;
			return (Pressure(rho + drho, u) - Pressure(rho - drho, u)) / (2.0 * drho);
		}
		inline type dPdu(const type rho, const type u) const{
			const type du = 0.0001;
			return (Pressure(rho, u + du) - Pressure(rho, u - du)) / (2.0 * du);
		}
		public:
		Tillotson(const type a_rho0, const type a_u0, const type a_uiv, const type a_ucv, const type a_A, const type a_B, const type a_a, const type a_b, const type a_alpha, const type a_beta){
            //in MKS unit...
            //From Brundage 2013
            //"Implementation of Tillotson Equation of State for Hypervelocity Impact of Metals, Geologic Materials, and Liquids"
            rho0  = a_rho0; // Density                          [kg/m^3]
            u0    = a_u0;   // Initial Energy                   [J/kg]
            uiv   = a_uiv;  // Energy at incipient vaporization [J/kg]
            ucv   = a_ucv;  // Energy at complete vaporization  [J/kg]
            A     = a_A;    // Bulk modulus                     [Pa]
            B     = a_B;    // Tillotson parameter              [Pa]
            a     = a_a;    // Tillotson parameter              [dimension-less]
            b     = a_b;    // Tillotson parameter
            alpha = a_alpha;// Tillotson parameter
            beta  = a_beta; // Tillotson parameter
		}
		inline type Pressure(const type dens, const type eng) const{
			const type p_min = 1.0e+7;
			if(dens >= rho0 || eng < uiv){
				return std::max(P_co(dens, eng), p_min);
			}else if(dens < rho0 && eng > ucv){
				return std::max(P_ex(dens, eng), p_min);
			}else{
				return std::max(((eng - uiv) * P_ex(dens, eng) + (ucv - eng) * P_co(dens, eng)) / (ucv - uiv), p_min);
			}
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			return sqrt(std::max(Pressure(dens, eng) / (dens * dens) * dPdu(dens, eng) + dPdrho(dens, eng), 0.0) + 1.0e-16);
		}
	};

	template <typename type> class ANEOS : public EoS_t<type>{
		/**
		 * A member variable that contains all the EoS data from ANEOS. The lines
		 * represent density, the columns energy. The data values for each entry
		 * are
		 * Density (kg/m3) Temperature (K) Energy (J/kg) Pressure (Pa) Sound speed (m/s) Entropy (J/kg/K).
		 */
		std::vector<std::vector<std::array<type, 6> > > eos_data;

		/**
		 * Linear vectors of the values of the corresponding fields in eos_data. Stored
		 * for making finding correct lines in eos_data simpler.
		 */
		std::vector<double> densities;
		std::vector<double> energies;

		/**
		 * A function that returns the corresponding table indices for a given
		 * density and energy. Note that the current implementation is very naive (but short).
		 * If the data
		 * points were equidistant we do not even need a search. But before optimizing we
		 * should first check if this actually matters in terms of performance.
		 * The function returns the table index of the next larger density and energy.
		 */
		std::pair<unsigned int, unsigned int>
		get_table_index(const type density, const type energy) const{
			auto line = std::lower_bound(densities.begin(),densities.end(),density);
			const unsigned int line_index = std::distance(densities.begin(),line);
			auto column = std::lower_bound(energies.begin(),energies.end(),energy);
			const unsigned int column_index = std::distance(energies.begin(),column);

			const unsigned int max_line_index = densities.size()-1;
			const unsigned int max_column_index = energies.size()-1;

			return std::make_pair(std::min(line_index,max_line_index),std::min(column_index,max_column_index));
		}

		const type
		get_interpolated_value (const type density,
				const type energy,
				const unsigned int property_index) const{
			const std::pair<unsigned int, unsigned int> data_index = get_table_index(density, energy);

			// If we are not at the boundaries of the data table
			if ((data_index.first < densities.size()-1) &&
					(data_index.second < energies.size()-1))
			{
				// compute the interpolation weights of this density and energy
				const type xi = (density - densities[data_index.first-1]) /
						(densities[data_index.first] - densities[data_index.first-1]);

				const type eta = (energy - energies[data_index.second-1]) /
						(energies[data_index.second] - energies[data_index.second-1]);

				// use these coordinates for a bilinear interpolation
				return  (1-xi)*(1-eta) * eos_data[data_index.first-1][data_index.second-1][property_index] +
						xi    *(1-eta) * eos_data[data_index.first]  [data_index.second-1][property_index] +
						(1-xi)*eta     * eos_data[data_index.first-1][data_index.second]  [property_index] +
						xi    *eta     * eos_data[data_index.first]  [data_index.second]  [property_index];
			}
			else
			{
				// Return the boundary value
				return eos_data[data_index.first][data_index.second][property_index];
			}
			return eos_data[data_index.first][data_index.second][property_index];
		}

		public:
		/**
		 * Construct the EoS by reading its data from a file. The file format looks
		 * as follows:
		 * TODO fill in format description, check that my implementation is correct
		 */
		ANEOS(const std::string &filename){
			std::ifstream input;
			input.open(filename, std::ios::in);
			if(input.fail() == true){
				std::cout << "Cannot open file..." << filename << std::endl;
				return;
			}

			const unsigned int n_expected_fields = 5;
			unsigned int n_expected_lines = 0;
			unsigned int n_expected_columns = 0;

			unsigned int n_lines = 0;
			unsigned int column_index = 0;
			std::string line;

			while(std::getline(input, line)){
				// read the table size, which has to be of the format # n_rows n_cols
				if (n_lines == 0 && n_expected_lines == 0 && n_expected_columns == 0){
					std::string tmp;
					std::istringstream stream(line);

					stream >> tmp;
					stream >> n_expected_lines;
					stream >> n_expected_columns;

					++n_lines;
					column_index = 0;
					eos_data.push_back(std::vector<std::array<type,6> >());
				}

				// skip other comments
				if(line[0] == '#')
					continue;

				// if we have reached the expected number of columns add a new row
				if(column_index == n_expected_columns){
					++n_lines;
					column_index = 0;
					eos_data.push_back(std::vector<std::array<type,6> >());
					}

				const unsigned int line_index = n_lines - 1;
				eos_data[line_index].push_back(std::array<type,6>());

				// make stream for reading
				std::istringstream stream(line);
				unsigned int field_index = 0;
				type tmp;

				while(stream >> tmp)
				{
					eos_data[line_index][column_index][field_index] = tmp;

					if (column_index == 0 && field_index == 0)
						densities.push_back(tmp);
					if (line_index == 0 && field_index == 1)
						energies.push_back(tmp);

					++field_index;
				}
				++column_index;
			}
			input.close();

			test_data();
		}

		inline type Pressure(const type dens, const type eng) const{
			return get_interpolated_value(dens,eng,3);
		}

		inline type SoundSpeed(const type dens, const type eng) const{
			return get_interpolated_value(dens,eng,4);
		}

		void test_data() const{
			std::ofstream output;
			output.open("test_output.txt", std::ios::out);
			output.setf(std::ios::scientific);
			output.precision(8);
			if(output.fail() == true){
				std::cout << "Cannot open file test_output.txt" << std::endl;
				return;
			}

			const unsigned int n_expected_fields = 5;
			unsigned int n_expected_lines = densities.size();
			unsigned int n_expected_columns = energies.size();

			for (unsigned int line = 0; line < n_expected_lines; ++line)
				for (unsigned int column = 0; column < n_expected_columns; ++column)
				{
					output << ' ';
					for (unsigned int field = 0; field < n_expected_fields; ++field)
					{
						output << eos_data[line][column][field] << "  ";
					}

					output << std::endl;
				}
			output.close();
		}
	};
}

static const EoS::IdealGas<PS::F64>  Monoatomic(5./3.);
static const EoS::IdealGas<PS::F64>  Diatomic  (1.4);
static const EoS::Tillotson<PS::F64> Granite   (2680.0, 16.0e+6, 3.5e+6, 18.00e+6,  18.0e+9,  18.0e+9, 0.5, 1.3, 5.0, 5.0);
static const EoS::Tillotson<PS::F64> Iron      (7800.0,  9.5e+6, 2.4e+6 , 8.67e+6, 128.0e+9, 105.0e+9, 0.5, 1.5, 5.0, 5.0);
static const EoS::ANEOS<PS::F64> AGranite      ("granite.rho_u.txt");

