#include <iostream>
#include <random>

int main()
{
	std::random_device r;
	std::default_random_engine e1(r());
	std::uniform_int_distribution<int> uniform_int(1, 20);
	std::uniform_real_distribution<double> uniform_real(1, 20);

	std::cout << "polyhedron([";

	auto num_points = uniform_int(e1);
	for(auto i=0; i<num_points; ++i)
	{
		if(i>0)
			std::cout << ',';
		std::cout << '[' << uniform_real(e1) << ',' << uniform_real(e1) << ',' << uniform_real(e1) << ']';
	}

	std::cout << "],[";

	std::uniform_int_distribution<int> uniform_face(3, 6);
	auto num_faces = uniform_int(e1);
	for(auto f=0; f<num_faces; ++f)
	{
		if(f>0)
			std::cout << ',';
		std::cout << '[';
		auto num_indices = uniform_face(e1);
		for(auto i=0; i<num_indices; ++i) {
			std::uniform_int_distribution<int> uniform_index(0, num_points);
			if(i>0)
				std::cout << ',';
			std::cout << uniform_index(e1);
		}
		std::cout << ']';
	}
	std::cout << "]);";
}
