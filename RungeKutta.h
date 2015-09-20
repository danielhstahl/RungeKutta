#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream> //debugging


class RungeKutta { //generic class, can take complex numbers, multivariate functions, etc
	private:
		int numSteps;
		double t;
	public:
		RungeKutta(double,  int);
		template<typename T, typename FN, typename... ARGS>
		std::vector<T> compute(FN&& fn, std::vector<T> initialValues, ARGS&&... args){
			double h=t/numSteps;
			double hlfh=h*.5;
			double sixthh=h/6.0;
			int n=initialValues.size();
			for(int i=0; i<numSteps; i++){
				std::vector<T> newVal(n);
				std::vector<T> k1=fn(h*i, initialValues, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k1[j]*hlfh;
				}
				std::vector<T> k2=fn(h*i+hlfh, newVal, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k2[j]*hlfh;
				}
				std::vector<T> k3=fn(h*i+hlfh, newVal, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k3[j]*h;
				}
				std::vector<T> k4=fn(h*i+h, newVal, args...);
				for(int j=0; j<n;j++){
					initialValues[j]=initialValues[j]+(k1[j]+k2[j]*2+k3[j]*2+k4[j])*sixthh;
				}
			}
			return initialValues;
		}
};
#endif
