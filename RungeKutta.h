#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream> //debugging
#include "FunctionalUtilities.h" //this is a personal repository

namespace rungekutta { //generic class, can take complex numbers etc
	/*private:
		int numSteps;
		double t;
		double h;
		double hlfh;
		double sixthh;
	public:
		RungeKutta(double t_, int numSteps_){
			t=t_;
			numSteps=numSteps_;
			h=t/numSteps;
			hlfh=h*.5;
			sixthh=h/6.0;
		}*/
	template<typename T, typename FN>
	std::vector<T> compute(const auto& t, const auto& numSteps, const std::vector<T>& initialValues, FN&& fn){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		auto sixthh=h/6.0;
		int n=initialValues.size();
		std::vector<T> newVal(initialValues); //temporarily hold values
		for(int i=0; i<numSteps; ++i){
			//std::vector<T> k1=fn(h*i, newVal);

			//initialValues

			std::vector<T> k1=fn(h*i, newVal);
			for(int j=0; j<n;++j){
				newVal[j]=newVal[j]+k1[j]*hlfh;
			}
			std::vector<T> k2=fn(h*i+hlfh, newVal);
			for(int j=0; j<n;++j){
				newVal[j]=newVal[j]+k2[j]*hlfh-k1[j]*hlfh;
			}
			std::vector<T> k3=fn(h*i+hlfh, newVal);
			for(int j=0; j<n;++j){
				newVal[j]=newVal[j]+k3[j]*h-k2[j]*hlfh;
			}
			std::vector<T> k4=fn(h*i+h, newVal);
			for(int j=0; j<n;++j){
				newVal[j]=newVal[j]+sixthh*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j]);
			}
		}
		return newVal;
	}

	/*auto functionalRG(auto&& fn){
		return [&](const auto& t, const auto& h, const auto& hlfh, auto&& y){
			return [&](auto&& dy1){
				return [&](auto&& dy2){
					return [&](auto&& dy3){
						return [&](auto&& dy4){
							return futilities::for_each_parallel(dy4, [&](const auto& val, const auto& index){
								return y[index]+(dy1[index]+2.0*dy2[index]+2.0*dy3[index]+val)/6.0;
							});
						}(fn(t+h, futilities::for_each_parallel(dy3, [&](const auto& val, const auto& index){return y[index]+val*h;})));
					}(fn(t+hlfh, futilities::for_each_parallel(dy2, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
				}(fn(t+hlfh, futilities::for_each_parallel(dy1, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
			}(fn(t, y)); //evaluates with this argument
		};
	}*/
	template<typename Number>
	std::vector<Number> computeFunctional(const auto& t, const auto& numSteps, const std::vector<Number>& initialValues, auto&& fn){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		//auto myResult=initialValues;
		auto fnc=[&](const auto& t, const auto& h, const auto& hlfh, const auto& y){
			return [&](const auto& dy1){
				return [&](const auto& dy2){
					return [&](const auto& dy3){
						return [&](const auto& dy4){
							return futilities::for_each_parallel_copy(dy4, [&](const auto& val, const auto& index){
								return y[index]+h*(dy1[index]+2.0*dy2[index]+2.0*dy3[index]+val)/6.0;
							});
						}(fn(t+h, futilities::for_each_parallel_copy(dy3, [&](const auto& val, const auto& index){return y[index]+val*h;})));
					}(fn(t+hlfh, futilities::for_each_parallel_copy(dy2, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
				}(fn(t+hlfh, futilities::for_each_parallel_copy(dy1, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
			}(fn(t, y)); //evaluates with this argument
		};
		return futilities::recurse(numSteps, initialValues, [&](const auto& val, const auto& index){
			return fnc(index*h, h, hlfh, val);
		});
	}
}
#endif
