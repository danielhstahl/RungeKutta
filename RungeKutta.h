#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream> //debugging


class RungeKutta { //generic class, can take complex numbers etc
	private:
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
		}
		template<typename T, typename FN>
		std::vector<T> compute(FN&& fn, const std::vector<T>& initialValues){
			int n=initialValues.size();
			std::vector<T> newVal(initialValues); //temporarily hold values
			for(int i=0; i<numSteps; i++){
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
					newVal[j]=newVal[j]+(k1[j]+k2[j]*2+k3[j]*2+k4[j])*sixthh;
				}
			}
			return newVal;
		}
};
#endif
