#ifndef CHEMMISOL_HOMOTOPY_H
#define CHEMMISOL_HOMOTOPY_H

#include "newton.h"

namespace chemmisol {

	template<typename X>
	class H {
		private:
			double t;
			const std::function<X(const X&)>& f; // target system
			const std::function<X(const X&)>& g; // start system
		public:
			H(
					double t,
					const std::function<X(const X&)>& f,
					const std::function<X(const X&)>& g)
				: t(t), f(f), g(g) {
				}

			X operator()(const X& x) {
				return typename X::value_type(1-t) * g(x)
					+ typename X::value_type(t) * f(x);
			}
	};

	template<typename X, typename M>
	class dH {
		private:
			double t;
			const std::function<M(const X&)>& df;
			const std::function<M(const X&)>& dg;
		public:
			dH(
					double t,
					const std::function<M(const X&)>& df,
					const std::function<M(const X&)>& dg)
				: t(t), df(df), dg(dg) {
				}

			M operator()(const X& x) {
				return typename X::value_type(1-t) * dg(x)
					+ typename X::value_type(t) * df(x);
			}
	};

	template<typename X, typename M>
		class Homotopy {
			private:
				std::list<X> x0;
				std::function<X(const X&)> f; // target system
				std::function<M(const X&)> df;
				std::function<X(const X&)> g; // start system
				std::function<M(const X&)> dg;

				X solve_iter(double t, double delta_t, const X& x, std::size_t n) const;

			public:
				Homotopy(
						const std::list<X>& x0,
						std::function<X(const X&)> f,
						std::function<M(const X&)> df,
						std::function<X(const X&)> g,
						std::function<M(const X&)> dg
						) :
					x0(x0), f(f), df(df), g(g), dg(dg) {
					}

				std::list<X> solve_iter(std::size_t n) const;
		};

	template<typename X, typename M>
		std::list<X> Homotopy<X, M>::solve_iter(std::size_t n) const {
			std::list<X> results;
			for(const auto& x : x0) {
				double delta = 1.0/n;
				results.push_back(solve_iter(delta, delta, x, n));
			}
			return results;
		}

	template<typename X, typename M>
		X Homotopy<X, M>::solve_iter(
				double t, double delta_t, const X& x, std::size_t n
				) const {
			CHEM_LOG(TRACE) << "[HOMOTOPY] Current X: " << x << " (t=" << t << ")";
			Newton<X, M> newton(
					x,
					H<X>(t, f, g),
					dH<X, M>(t, df, dg));
			
			X xt = newton.solve_iter(n);
			if(t < 1.0) {
				t = t + delta_t;
				return solve_iter(t, delta_t, xt, n);
			}
			return xt;
		}
}
#endif
