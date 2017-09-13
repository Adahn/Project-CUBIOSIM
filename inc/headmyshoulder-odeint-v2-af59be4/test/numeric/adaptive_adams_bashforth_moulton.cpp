/* Boost numeric test of the adams-bashforth steppers test file

 Copyright 2013 Karsten Ahnert
 Copyright 2013-2015 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

// disable checked iterator warning for msvc
#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE numeric_adaptive_adams_bashforth_moulton

#include <iostream>
#include <cmath>

#include <boost/array.hpp>

#include <boost/test/unit_test.hpp>

#include <boost/mpl/vector.hpp>

#include <boost/numeric/odeint.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;
namespace mpl = boost::mpl;

typedef double value_type;

typedef boost::array< double , 2 > state_type;
typedef runge_kutta_fehlberg78<state_type> initializing_stepper;

// harmonic oscillator, analytic solution x[0] = sin( t )
struct osc
{
    void operator()( const state_type &x , state_type &dxdt , const double t ) const
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
    }
};

BOOST_AUTO_TEST_SUITE( numeric_adaptive_adams_bashforth_moulton_test )


/* generic test for all adams bashforth steppers */
template< class Stepper >
struct perform_adaptive_adams_bashforth_moulton_test
{
    void operator()( void )
    {
        Stepper stepper;
        initializing_stepper init_stepper;

        const int o = stepper.order()+1; //order of the error is order of approximation + 1

        const state_type x0 = {{ 0.0 , 1.0 }};
        state_type x1 = x0;
        double t = 0.0;
        double dt = 0.25;
        // initialization, does a number of steps to self-start the stepper with a small stepsize
        stepper.initialize( init_stepper, osc() , x1 , t ,  dt);
        double A = std::sqrt( x1[0]*x1[0] + x1[1]*x1[1] );
        double phi = std::asin(x1[0]/A) - t;
        // more steps necessary to "counteract" the effect from the lower order steps
        for( size_t n=0 ; n < stepper.steps ; ++n )
        {
            stepper.do_step( osc() , x1 , t , dt );
            t += dt;
        }
        // now we do the actual step
        stepper.do_step( osc() , x1 , t , dt );
        // only examine the error of the adams-bashforth step, not the initialization
        const double f = 2.0 * std::abs( A*sin(t+dt+phi) - x1[0] ) / std::pow( dt , o ); // upper bound
        
        std::cout << o << " , " << f << std::endl;

        /* as long as we have errors above machine precision */
        while( f*std::pow( dt , o ) > 1E-16 )
        {
            x1 = x0;
            t = 0.0;
            stepper.initialize( init_stepper, osc() , x1 , t , dt );
            A = std::sqrt( x1[0]*x1[0] + x1[1]*x1[1] );
            phi = std::asin(x1[0]/A) - t;
            // now we do the actual step
            stepper.do_step( osc() , x1 , t , dt );
            stepper.reset();
            // only examine the error of the adams-bashforth step, not the initialization
            std::cout << "Testing dt=" << dt << " , " << std::abs( A*sin(t+dt+phi) - x1[0] ) << std::endl;
            BOOST_CHECK_LT( std::abs( A*sin(t+dt+phi) - x1[0] ) , f*std::pow( dt , o ) );
            dt *= 0.5;
        }
    }
};

typedef mpl::vector<
    adaptive_adams_bashforth_moulton< 2 , state_type > ,
    adaptive_adams_bashforth_moulton< 3 , state_type > ,
    adaptive_adams_bashforth_moulton< 4 , state_type > ,
    adaptive_adams_bashforth_moulton< 5 , state_type > ,
    adaptive_adams_bashforth_moulton< 6 , state_type > ,
    adaptive_adams_bashforth_moulton< 7 , state_type >
    > adaptive_adams_bashforth_moulton_steppers;

BOOST_AUTO_TEST_CASE_TEMPLATE( adaptive_adams_bashforth_moulton_test , Stepper, adaptive_adams_bashforth_moulton_steppers )
{
    perform_adaptive_adams_bashforth_moulton_test< Stepper > tester;
    tester();
}

BOOST_AUTO_TEST_SUITE_END()
