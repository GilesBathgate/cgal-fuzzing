#define BOOST_DISABLE_ASSERTS
#include <iostream>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#include <CGAL/Modifier_base.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <cmath>
#include <vector>
#include <tuple>
#include <unordered_set>
#include <set>

using Point = std::vector<double>;
using Points = std::vector<Point>;
using Facet = std::vector<Points::size_type>;
using Facets = std::vector<Facet>;
using Polyhedron = std::tuple<Points,Facets>;
using Components = std::vector<Polyhedron>;

using K=CGAL::Exact_predicates_exact_constructions_kernel;
using HalfedgeDS=CGAL::Polyhedron_3<K>::HalfedgeDS;

/*
static bool valid_indexing(const Points& points,const Facets& facets)
{
	std::unordered_set<int> used;
	std::set<CGAL::Point_3<K>> unique;
	for(const auto& face: facets) {
		for(const auto& index: face) {
			const auto& p=points[index];
			unique.insert(CGAL::Point_3<K>(p[0],p[1],p[2]));
			used.insert(index);
		}
	}

	// used == unique == points
	return used.size()==unique.size()&&used.size()==points.size();
}
*/

class Builder : public CGAL::Modifier_base<HalfedgeDS>
{
	const Polyhedron& poly;
	bool complete;

public:
	Builder(const Polyhedron& p) : poly(p), complete(false) {}

	void operator()(HalfedgeDS& hds) override {

		const auto& points = std::get<Points>(poly);
		const auto& facets = std::get<Facets>(poly);

		//sanity checks
		for(const auto& point: points) {
			if(point.size()!=3) return;
			for(const auto& v: point) {
				if(!std::isfinite(v)) return;
			}
		}

		for(const auto& face: facets) {
			for(const auto& index: face) {
				if(index<0||index>=points.size())
					return;
			}
		}

		CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds,true);
		builder.begin_surface(points.size(),facets.size());
		for(const auto& p: points)
			builder.add_vertex(CGAL::Point_3<K>(p[0],p[1],p[2]));

		for(const auto& facet: facets) {
			if(!builder.test_facet(facet.begin(),facet.end())) {
				builder.rollback();
				return;
			}
			builder.add_facet(facet.begin(),facet.end());
		}

		builder.end_surface();
		builder.remove_unconnected_vertices();

		for(auto f=hds.faces_begin(); f!=hds.faces_end(); ++f) {
			if(f->facet_degree()>3) continue;
			auto h1=f->halfedge(),h2=h1->next(),h3=h2->next();
			if(CGAL::collinear(h1->vertex()->point(),h2->vertex()->point(),h3->vertex()->point()))
				return;
		}

		//check for zero length edges
		for(auto h=hds.halfedges_begin(); h!=hds.halfedges_end(); ++h) {
			if(0.0 == CGAL::squared_distance(h->vertex()->point(),h->opposite()->vertex()->point()))
				return;
		}

		for(auto h=hds.faces_begin(); h!=hds.faces_end(); ++h) {
			if(h->facet_degree()>3) {
				//check for self-intersecting facets
				auto b=h->facet_begin(),e(b),o(b);
				CGAL_For_all(b,e) {
					CGAL_For_all(o,e) {
						if(b==o||b->next()==o||o->next()==b||b->prev()==o||o->prev()==b) continue;
						CGAL::Segment_3<K> sb(b->vertex()->point(),b->opposite()->vertex()->point());
						CGAL::Segment_3<K> so(o->vertex()->point(),o->opposite()->vertex()->point());
						if(CGAL::do_intersect(sb,so)) return;
					}
				}

				//check for non-coplanar facets
				CGAL::Point_3<K> p[4];
				for(auto i=0; i<2; ++i,++b)
					p[i]=b->vertex()->point();
				CGAL_For_all(b,e) {
					p[2]=b->vertex()->point();
					if(!CGAL::collinear(p[0],p[1],p[2])) break;
				}
				if(b==e) return;
				e=b;
				CGAL_For_all(b,e) {
					p[3]=b->vertex()->point();
					if(!CGAL::coplanar(p[0],p[1],p[2],p[3])) return;
				}
			}
		}

		//check for facet fans with more than 2 border edges
		for(auto v=hds.vertices_begin(); v!=hds.vertices_end(); ++v) {
			if(v->vertex_degree()<=2) continue;
			unsigned borderEdges=0;
			auto b=v->vertex_begin(),e(b);
			CGAL_For_all(e,b) {
				if(e->is_border_edge())
					++borderEdges;
			}
			if(borderEdges>2)
				return;
		}

		complete = true;
	}

	bool is_complete() { return complete; }
};

using namespace boost::spirit::qi;

template <typename T>
using Rule = rule<std::string::const_iterator, T, space_type>;

void test(const std::string& input)
{
	Components result;

	Rule<Point> point_rule = '[' >> (double_ % ',') >> ']';
	Rule<Points> points_rule = '[' >> point_rule % ',' >> ']';
	Rule<Facet> facet_rule = '[' >> (int_ % ',') >> ']';
	Rule<Facets> facets_rule = '[' >> facet_rule % ',' >> ']';
	Rule<Polyhedron> polyhedron_rule = "polyhedron(" >> points_rule >> ',' >> facets_rule >> ");";
	Rule<Components> components_rule = "difference(){" >> *polyhedron_rule >> "}";

	bool parsed = phrase_parse(
			input.begin(),input.end(),
			components_rule,
			space,result);

	if(parsed && result.size()==2) {
		Builder b1(result[0]);

		CGAL::Polyhedron_3<K> P1;
		P1.delegate(b1);
		if(!b1.is_complete()) return;

		Builder b2(result[1]);

		CGAL::Polyhedron_3<K> P2;
		P2.delegate(b2);
		if(!b2.is_complete()) return;

		CGAL::Nef_polyhedron_3<K> N1(P1);
		CGAL::Nef_polyhedron_3<K> N2(P2);

		CGAL::Nef_polyhedron_3<K> N3 = N1.difference(N2);

		std::cout << "valid!" << std::endl;
	}
}


int main() {
	std::string input;
	std::getline(std::cin, input);
	test(input);
}

/*
#pragma clang optimize off
__AFL_FUZZ_INIT();
int main() {
#ifdef __AFL_HAVE_MANUAL_CONTROL
	__AFL_INIT();
#endif
	unsigned char* buf = __AFL_FUZZ_TESTCASE_BUF;
	while (__AFL_LOOP(10000)) {
		int len = __AFL_FUZZ_TESTCASE_LEN;
		if (len < 10) continue;
		test(std::string((char*)buf, len));
	}
	return 0;
}
*/
