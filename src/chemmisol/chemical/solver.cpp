#include "chemmisol/chemical/solver.h"
#include "chemmisol/math/newton.h"
#include "chemmisol/math/arithmetics.h"
#include <cassert>

namespace chemmisol {
	namespace solver {
		AbsoluteNewton default_solver;

		X AbsoluteNewton::solve(const ChemicalSystem& system) const {
			ReducedChemicalSystem<X> reduced_system(system);
			F<X, M> f(reduced_system, system);
			// Initial activities in the system
			X reduced_activities = reduced_system.reducedActivities();
			reduced_activities = chemmisol::AbsoluteNewton<X, M>(
					reduced_activities,
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(system.getMaxIteration());
			return reduced_system.completeActivities(reduced_activities);
		}

		X G::compute_degrees(
				const ReducedChemicalSystem<CX>& reduced_system,
				const ChemicalSystem& system
				) {
			X deg(reduced_system.xSize());
			for(std::size_t i = 0; i < reduced_system.reactionOffset(); i++)
				// Mass conservation law polynoms
				deg[i] = 1.0;
			for(const auto& reaction : system.getReactions())
				deg[reduced_system.reactionOffset()
					+ reaction->getIndex()] = system.degree(*reaction);
			return deg;
		}

		CX G::g(const CX& reduced_activities) const {
			CX gx(reduced_system.xSize());
			for(std::size_t i = 0; i < reduced_system.xSize(); i++) {
				gx[i] = a[i] * std::pow(reduced_activities[i], degrees[i]) - b[i];
			}
			return gx;
		}

		CM G::dg(const CX& reduced_activities) const {
			CM dg(reduced_system.xSize());
			for(std::size_t i = 0; i < reduced_system.xSize(); i++) {
				dg[i].resize(reduced_system.fxSize());
				// All other derivatives are equal to 0.0
				dg[i][i] = a[i] * degrees[i]
					* std::pow(reduced_activities[i], degrees[i]-1);
			}
			return dg;
		}

		std::list<CX> build_roots(
				const std::list<CX>& current_roots,
				const std::vector<std::vector<typename CX::value_type>>& roots,
				std::size_t i) {
			std::list<CX> local_roots;
			for(const CX& base : current_roots) {
				for(std::size_t j = 0; j < roots[i].size(); j++) {
					CX new_root(base);
					new_root[i] = roots[i][j];
					local_roots.push_back(new_root);
				}
			}
			if(i == roots.size()-1)
				return local_roots;
			return build_roots(local_roots, roots, i+1);
		}

		std::list<CX> build_roots(
				std::size_t n,
				const std::vector<std::vector<typename CX::value_type>>& roots) {
			std::list<CX> base;
			base.push_back(CX(n));
			return build_roots(
					base, roots, 0);
		}

		std::list<CX> G::initValues() const {
			std::vector<std::vector<typename CX::value_type>> _roots(reduced_system.xSize());
			for(std::size_t i = 0; i < reduced_system.xSize(); i++) {
				auto roots_list = roots(b[i] / a[i], degrees[i]);
				assert(roots_list.size() > 0);
				_roots[i] = {roots_list.begin(), roots_list.end()};
				CHEM_LOGV(8) << "roots[" << i << "]: " << roots_list;
			}
			return build_roots(reduced_system.xSize(), _roots);
		}
	}
}
