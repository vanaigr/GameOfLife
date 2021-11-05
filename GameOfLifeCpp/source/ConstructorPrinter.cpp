#include<iostream>

template<typename D>
class Printer {
public:
	Printer(std::string&& msg, D&& data) {
		std::cout << msg << data;
	}
};