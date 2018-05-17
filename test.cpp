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
TEST_CASE("Test functionalRG 2d", "[RG]"){
	double t=2.0;
	int numSteps=1024;
	std::vector<double> initialValues={1.0, 1.0};

	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::compute_efficient_2d(t, numSteps, std::move(initialValues), [](const auto& t, const auto& y1, const auto& y2){
		return std::vector<double>({y1*t, y2*t});
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test functionalRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues[0]==Approx(exp(2.0)));
}

TEST_CASE("Test functionalRG single", "[RG]"){
	double t=2.0;
	int numSteps=1024;
	double initialValues=1.0;
	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::computeFunctional(t, numSteps, initialValues, [](const auto& t, const auto& y){
		return y*t;
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test functionalRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues==Approx(exp(2.0)));
}
TEST_CASE("Test functionalRG move", "[RG]"){
	double t=2.0;
	int numSteps=1024;
	std::vector<double> initialValues={1.0, 1.0};

	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::computeFunctional_move(t, numSteps, std::move(initialValues), [](const auto& t, const auto& y){
		return std::vector<double>({y[0]*t, y[1]*t});
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test functionalRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues[0]==Approx(exp(2.0)));
}

TEST_CASE("Test functionalRG single move", "[RG]"){
	double t=2.0;
	int numSteps=1024;
	double initialValues=1.0;
	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::computeFunctional_move(t, numSteps, std::move(initialValues), [](const auto& t, const auto& y){
		return y*t;
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test functionalRG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	REQUIRE(initialValues==Approx(exp(2.0)));
}
TEST_CASE("Test speed of convergence", "[RG]"){
	double t=2.0;
	int numSteps=32; //pretty low

	std::vector<double> initialValues={1.0, 1.0};

	auto t1 = std::chrono::high_resolution_clock::now();
	initialValues=rungekutta::computeFunctional_move(t, numSteps, std::move(initialValues), [](const auto& t, const auto& y){
		return std::vector<double>({y[0]*t, y[1]*t});
	});   
	auto t2 = std::chrono::high_resolution_clock::now();
	REQUIRE(initialValues[0]==Approx(exp(2.0))); 
}

