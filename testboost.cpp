
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

int main()
{
    ublas::vector<double> v1(2), v2(4);
    for (unsigned i = 0; i < 2; ++i)
        v1 (i) = i;
     for (unsigned i = 0; i < 4; ++i)
        v2 (i) = i;   

    ublas::matrix<double> m(v1.size(),v1.size()*v2.size());
    m =  outer_prod(v1, v2);

    std::cout << v1 << '\n'
              << v2 << '\n'
              << m << '\n';

}