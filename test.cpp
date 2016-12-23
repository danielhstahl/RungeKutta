#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "RungeKutta.h"
#include <vector>
#include <chrono>

TEST_CASE("Test functionalRG", "[RG]"){
	double t=2.0;
	int numSteps=1024;
	std::vector<double> initialValues={1.0, 1.0};

	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::computeFunctional(t, numSteps, initialValues, [](const auto& t, const auto& y){
		return std::vector<double>({y[0]*t, y[1]*t});
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test functionalRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues[0]==Approx(exp(2.0)));
}
/*TEST_CASE("Test standard RG", "[RG]"){
	double t=2.0;
	int numSteps=1024; 
	std::vector<double> initialValues={1.0, 1.0};
	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::compute(t, numSteps, initialValues, [](const auto& t, const auto& y){  
		return std::vector<double>({y[0]*t, y[1]*t});
	});    
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test standardRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues[0]==Approx(exp(2.0)));
}
*/