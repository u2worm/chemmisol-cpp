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

				SolverResult<X> solve(
						const X& x, std::size_t i,
						std::size_t homotopy_n, std::size_t local_solver_n) const;

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

				std::list<SolverResult<X>> solve(std::size_t homotopy_n, std::size_t local_solver_n) const;
		};

	template<typename X, typename M>
		std::list<SolverResult<X>> Homotopy<X, M>::solve(
				std::size_t homotopy_n, std::size_t local_solver_n) const {
			std::list<SolverResult<X>> results;
			CHEM_LOG(INFO) << "[HOMOTOPY] Start exhaustive homotopy from " << x0.size() << " starting points.";
			std::size_t i = 0;
			for(const auto& x : x0) {
				CHEM_LOG(INFO) << "[HOMOTOPY]   i=" << i;
				CHEM_LOG(TRACE) << "[HOMOTOPY] Current X: " << x << " (t=0)";
				results.push_back(solve(x, 0, homotopy_n, local_solver_n));
				++i;
			}
			return results;
		}

	template<typename X, typename M>
		SolverResult<X> Homotopy<X, M>::solve(
				const X& x, std::size_t i, std::size_t homotopy_n, std::size_t local_solver_n
				) const {
			double t = ((double) i+1)/homotopy_n;
			Newton<X, M> newton(
					x,
					H<X>(t, f, g),
					dH<X, M>(t, df, dg));
			
			auto result = newton.solve_iter(local_solver_n);
			CHEM_LOG(TRACE) << "[HOMOTOPY] Current X: " << result.x << " (t=" << t << ")";
			if(i == homotopy_n-1 || !result.isFinite()) {
				CHEM_LOG(INFO) << "[HOMOTOPY] Final X: " << result.x << " f(X)=" << result.f_x;
				return result;
			}
			return solve(result.x, i+1, homotopy_n, local_solver_n);
		}
}
#endif
