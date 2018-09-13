#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

class ParamFile{
	std::map<std::string, std::string> param;
	public:
	ParamFile(const char* const file){
		std::ifstream input;
		input.open(file, std::ios::in);
		if(input.fail() == true){
			std::cout << "Cannot open file..." << std::endl;
			return;
		}
		std::string line;
		while(std::getline(input, line)){
			//comment out
			if(line[0] == '#') continue;
			//empty line
			if(line.empty() == true) continue;
			//make stream
			std::istringstream stream(line);
			std::string tmp;
			std::vector<std::string> keyval;
			while(std::getline(stream, tmp, '=')){
				keyval.push_back(tmp);
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
	std::string getValueOf(std::string key){
		return param[key];
	}
};
