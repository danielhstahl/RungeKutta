#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream> //debugging
#include "FunctionalUtilities.h" //this is a personal repository
#include <type_traits>

#ifndef __IS_VECTOR
#define __IS_VECTOR
template<typename T> struct is_vector : public std::false_type {};

template<typename T, typename A> struct is_vector<std::vector<T, A>> : public std::true_type {};
#endif

namespace rungekutta { //generic class, can take complex numbers etc

	template<typename Number, typename Number1, typename Index, typename FN>
	std::vector<Number> computeFunctional_move(const Number1& t, const Index& numSteps, std::vector<Number>&& initialValues, FN&& fn, std::true_type){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		auto sixthh=h/6.0;
		auto fnc=[&](const auto& t, const auto& h, const auto& hlfh, const auto& y){
			return [&](auto&& dy1){
				return [&](auto&& dy2){
					return [&](auto&& dy3){
						return [&](auto&& dy4){
							return futilities::for_each_parallel(std::move(dy4), [&](const auto& val, const auto& index){
								return y[index]+(dy1[index]+2.0*dy2[index]+2.0*dy3[index]+val)*sixthh;
							});
						}(fn(t+h, futilities::for_each_parallel_copy(dy3, [&](const auto& val, const auto& index){return y[index]+val*h;})));
					}(fn(t+hlfh, futilities::for_each_parallel_copy(dy2, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
				}(fn(t+hlfh, futilities::for_each_parallel_copy(dy1, [&](const auto& val, const auto& index){return y[index]+val*hlfh;})));
			}(fn(t, y)); //evaluates with this argument
		};
		return futilities::recurse_move(numSteps, std::move(initialValues), [&](const auto& val, const auto& index){
			return fnc(index*h, h, hlfh, val);
		});
	}


	template<typename Number, typename Number1, typename Index, typename FN>
	Number computeFunctional_move(const Number1& t, const Index& numSteps,  Number&& initialValues, FN&& fn, std::false_type){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		auto sixthh=h/6.0;
		auto fnc=[&](const auto& t, const auto& h, const auto& hlfh, const auto& y){
			return [&](const auto& dy1){
				return [&](const auto& dy2){
					return [&](const auto& dy3){
						return [&](const auto& dy4){
							return y+(dy1+2.0*dy2+2.0*dy3+dy4)*sixthh;
						}(fn(t+h, y+dy3*h));
					}(fn(t+hlfh, y+dy2*hlfh));
				}(fn(t+hlfh, y+dy1*hlfh));
			}(fn(t, y)); //evaluates with this argument
		};
		return futilities::recurse_move(numSteps, std::move(initialValues), [&](const auto& val, const auto& index){
			return fnc(index*h, h, hlfh, val);
		});
	}

	template<typename Number, typename Number1, typename Index, typename FN>
	Number computeFunctional_move(const Number1& t, const Index& numSteps,  Number&& initialValues, FN&& fn){
		return computeFunctional_move(t, numSteps, std::move(initialValues), fn, is_vector<Number>{});
	}
	template<typename Number, typename Number1, typename Index, typename FN>
	std::vector<Number> computeFunctional(const Number1& t, const Index& numSteps, const std::vector<Number>& initialValues, FN&& fn, std::true_type){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		auto sixthh=h/6.0;
		auto fnc=[&](const auto& t, const auto& h, const auto& hlfh, const auto& y){
			return [&](auto&& dy1){
				return [&](auto&& dy2){
					return [&](auto&& dy3){
						return [&](auto&& dy4){
							return futilities::for_each_parallel(std::move(dy4), [&](const auto& val, const auto& index){
								return y[index]+(dy1[index]+2.0*dy2[index]+2.0*dy3[index]+val)*sixthh;
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


	template<typename Number, typename Number1, typename Index, typename FN>
	Number computeFunctional(const Number1& t, const Index& numSteps,  const Number& initialValues, FN&& fn, std::false_type){
		auto h=t/numSteps;
		auto hlfh=h*.5;
		auto sixthh=h/6.0;
		auto fnc=[&](const auto& t, const auto& h, const auto& hlfh, const auto& y){
			return [&](const auto& dy1){
				return [&](const auto& dy2){
					return [&](const auto& dy3){
						return [&](const auto& dy4){
							return y+(dy1+2.0*dy2+2.0*dy3+dy4)*sixthh;
						}(fn(t+h, y+dy3*h));
					}(fn(t+hlfh, y+dy2*hlfh));
				}(fn(t+hlfh, y+dy1*hlfh));
			}(fn(t, y)); //evaluates with this argument
		};
		return futilities::recurse(numSteps, initialValues, [&](const auto& val, const auto& index){
			return fnc(index*h, h, hlfh, val);
		});
	}

	template<typename Number, typename Number1, typename Index, typename FN>
	Number computeFunctional(const Number1& t, const Index& numSteps,  const Number& initialValues, FN&& fn){
		return computeFunctional(t, numSteps, initialValues, fn, is_vector<Number>{});
	}
}
#endif
