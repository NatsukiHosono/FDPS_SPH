#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <stdexcept>

/**
 * This is a utility function based on an identical function in the deal.II project
 * (https://github.com/dealii/dealii/blob/master/source/base/utilities.cc). It removes
 * whitespace at the beginning and end of a string.
 */
std::string
trim(const std::string &input)
{
std::string::size_type left = 0;
std::string::size_type right = input.size() > 0 ? input.size() - 1 : 0;

for (; left < input.size(); ++left)
  {
	if (!std::isspace(input[left]))
	  {
		break;
	  }
  }

for (; right >= left; --right)
  {
	if (!std::isspace(input[right]))
	  {
		break;
	  }
  }

return std::string(input, left, right - left + 1);
}

/**
 * This class handles reading in model input parameters from input files.
 */
class ParameterFile{
	std::map<std::string, std::string> param;
	public:
	ParameterFile(const std::string &file){
		std::ifstream input;
		input.open(file, std::ios::in);
		if(input.fail() == true){
			std::cout << "Cannot open file, assuming default values ..." << std::endl;
			return;
		}
		std::string line;
		while(std::getline(input, line)){
			// Remove comments
			auto trailing_comment = std::find(line.begin(),line.end(),'#');
			if (trailing_comment != line.end()){
				line.erase(trailing_comment,line.end());
			}
			line = trim(line);
			// Skip empty lines
			if(line.empty() == true) continue;
			// Parse line
			std::istringstream stream(line);
			std::string tmp;
			std::vector<std::string> keyval;
			while(std::getline(stream, tmp, '=')){
				keyval.push_back(trim(tmp));
			}
			param.insert(std::map<std::string, std::string>::value_type(keyval[0], keyval[1]));
		}
		input.close();
	}
	void show(){
		for(std::map<std::string, std::string>::iterator it = param.begin() ; it != param.end() ; ++ it){
			std::cout << it->first << " = " << it->second << std::endl;
		}
	}
	template <typename T>
	T getValueOf(const std::string &key, const T &default_parameter_value){

		T parameter_value = default_parameter_value;

		if (param.find(key) != param.end()){
			parameter_value = T();
			std::stringstream convert_from_string_to_value(param[key]);
			convert_from_string_to_value >> parameter_value;

		    if(convert_from_string_to_value.fail())
		    {
		      std::cout << "Convert error: cannot convert string '" << param[key] << "' to value" << std::endl;
		      throw;
		    }
		}
		else{
			std::cout << "Parameter <" << key << "> not found and no default value provided. Aborting." << std::endl;
			throw;
		}

		return parameter_value;
	  }
};
