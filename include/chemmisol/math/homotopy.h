#ifndef CHEMMISOL_HOMOTOPY_H
#define CHEMMISOL_HOMOTOPY_H

#include <omp.h>
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
				std::vector<X> x0;
				std::function<X(const X&)> f; // target system
				std::function<M(const X&)> df;
				std::function<X(const X&)> g; // start system
				std::function<M(const X&)> dg;

				SolverResult<X> solve(
						const X& x, std::size_t i,
						std::size_t homotopy_n, std::size_t local_solver_n) const;

			public:
				Homotopy(
						const std::vector<X>& x0,
						std::function<X(const X&)> f,
						std::function<M(const X&)> df,
						std::function<X(const X&)> g,
						std::function<M(const X&)> dg
						) :
					x0(x0), f(f), df(df), g(g), dg(dg) {
					}

				std::vector<SolverResult<X>> solve(std::size_t homotopy_n, std::size_t local_solver_n) const;
		};

	template<typename X, typename M>
		std::vector<SolverResult<X>> Homotopy<X, M>::solve(
				std::size_t homotopy_n, std::size_t local_solver_n) const {
			std::vector<SolverResult<X>> results(x0.size());
			CHEM_LOG(INFO) << "[HOMOTOPY] Start exhaustive homotopy from " << x0.size() << " starting points.";
#ifdef CHEMMISOL_OPENMP
			std::size_t n = 0;
#endif
			#pragma omp parallel for shared(results, x0, n)
			for(std::size_t i = 0; i < x0.size(); i++) {
#ifdef CHEMMISOL_OPENMP
				std::size_t _n;
				#pragma omp critical
				_n = ++n;

				CHEM_LOG(INFO) << "[HOMOTOPY]   i=" << i << " (thread " << omp_get_thread_num() << ", " << _n << "/" << x0.size() << ")";
#else
				CHEM_LOG(INFO) << "[HOMOTOPY]   i=" << i;
#endif
				CHEM_LOG(TRACE) << "[HOMOTOPY] Current X: " << x0[i] << " (t=0)";

				auto result =  solve(x0[i], 0, homotopy_n, local_solver_n);
				results[i] = result;
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
