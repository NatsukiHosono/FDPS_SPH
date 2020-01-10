#pragma once

#include <sstream>
#include "interpolation.h"
//#include "io.h"

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
		virtual type InternalEnergy(const type dens, const type ent) const = 0;
        virtual type Entropy(const type dens, const type eng) const = 0;
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

        inline type InternalEnergy(const type dens, const type ent) const{
            return 0;
        }

        inline type Entropy(const type dens, const type eng) const{
            return 0;
        }

		inline type HeatCapacityRatio() const{
			return hcr;
		}
	};

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

        inline type InternalEnergy(const type dens, const type ent) const{
            return 0;
        }

        inline type Entropy(const type dens, const type eng) const{
            return 0;
        }
	};

	template <typename type> class ANEOS : public EoS_t<type>{

		/**
		 * Linear vectors of the values of the corresponding fields in eos_data. Stored
		 * for making finding correct lines in eos_data simpler.
		 */
		std::vector<double> densities;
		std::vector<double> energies;
        std::vector<double> entropies;

        std::vector<double> full_densities;
        std::vector<double> full_energies;
        std::vector<double> full_entropies;

		public:
		/**
		 * Construct the EoS by reading its data from a file. The file format looks
		 * as follows:
		 * TODO fill in format description, check that my implementation is correct
		 */
        std::vector<std::vector<std::array<type, 6> > > eos_data;

		ANEOS(const std::string &filename){

            class readANEOSfile {
            public:
                static std::vector<std::vector<std::array<double, 6>>> readfile(
                        std::vector<double> &var1_vector,
                        std::vector<double> &var2_vector,
                        std::vector<double> &var3_vector,
                        std::vector<double> &var1_vector_full,
                        std::vector<double> &var2_vector_full,
                        std::vector<double> &var3_vector_full,
                        const std::string &file_path,
                        const int val1_property_index,
                        const int val2_property_index,
                        const int val3_property_index
                ) {

                    std::vector<std::vector<std::array<double, 6>>> eos_data;

                    const unsigned int n_expected_fields = 5;
                    unsigned int n_expected_lines = 0;
                    unsigned int n_expected_columns = 0;
                    unsigned int n_lines = 0;
                    unsigned int column_index = 0;

                    std::string();
                    std::ifstream input;
                    input.open(file_path, std::ios::in);

                    if (input.is_open()) {
                        std::string line;
                        while (std::getline(input, line)) {
                            // read the table size, which has to be of the format # n_rows n_cols
                            if (n_lines == 0 && n_expected_lines == 0 && n_expected_columns == 0) {
                                std::string tmp;
                                std::istringstream stream(line);

                                stream >> tmp;
                                stream >> n_expected_lines;
                                stream >> n_expected_columns;

                                ++n_lines;
                                column_index = 0;
                                eos_data.emplace_back();
                            }

                            // skip other comments
                            if (line[0] == '#')
                                continue;

                            // if we have reached the expected number of columns add a new row
                            if (column_index == n_expected_columns) {
                                ++n_lines;
                                column_index = 0;
                                eos_data.emplace_back();
                            }

                            const unsigned int line_index = n_lines - 1;
                            eos_data[line_index].push_back(std::array<double, 6>());

                            // make stream for reading
                            std::istringstream stream(line);
                            unsigned int field_index = 0;
                            double tmp;

                            while (stream >> tmp) {
                                eos_data[line_index][column_index][field_index] = tmp;

                                if (column_index == 0 && field_index == val1_property_index)
                                    var1_vector.push_back(tmp);
                                if (line_index == 0 && field_index == val2_property_index)
                                    var2_vector.push_back(tmp);
                                if (line_index == 0 && field_index == val3_property_index)
                                    var3_vector.push_back(tmp);

                                if (field_index == val1_property_index)
                                    var1_vector_full.push_back(tmp);
                                if (field_index == val2_property_index)
                                    var2_vector_full.push_back(tmp);
                                if (field_index == val3_property_index)
                                    var3_vector_full.push_back(tmp);

                                ++field_index;
                            }
                            ++column_index;
                        }
                    }
                    input.close();
                    return eos_data;
                }
            };

            eos_data = readANEOSfile::readfile(densities, energies,
                    entropies, full_densities, full_energies, full_entropies, filename, 0, 1, 5);
		}

		inline type Pressure(const type dens, const type eng) const{
			return BilinearInterpolation::interpolate(dens, eng, densities, energies, 3, eos_data);
		}

		inline type SoundSpeed(const type dens, const type eng) const{
            return BilinearInterpolation::interpolate(dens, eng, densities, energies, 4, eos_data);
		}

        inline type InternalEnergy(const type dens, const type ent) const{
            return EnergyInterpolation::interpolate(dens, ent, full_densities, full_entropies, full_energies, 1, eos_data, 120);
        }

        inline type Entropy(const type dens, const type eng) const{
            return BilinearInterpolation::interpolate(dens, eng, densities, energies, 5, eos_data);
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
static const EoS::ANEOS<PS::F64> AGranite      ("eos/granite.rho_u.txt");

