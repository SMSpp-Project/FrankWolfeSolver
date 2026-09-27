/*--------------------------------------------------------------------------*/
/*---------------------------- File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of FrankWolfeSolver on Blocks of the core alone.
 *
 * The father is an AbstractBlock whose Objective is a separable quadratic
 * DQuadFunction over the ColVariable of its K sub-Block; each sub-Block is an
 * AbstractBlock with a BoxConstraint on each of its ColVariable, a linear
 * Objective, and a BoxSolver as its oracle. The whole problem is then
 *
 *     min sum_i a_i x_i^2 + ( b_i + c_i ) x_i   s.t.   l_i <= x_i <= u_i
 *
 * whose optimum is known in closed form, x_i being the unconstrained minimum
 * - ( b_i + c_i ) / ( 2 a_i ) projected on [ l_i , u_i ]. Every algorithm of
 * FrankWolfeSolver (vanilla, Away-step, Blended Pairwise), under both of its
 * bookkeeping modes (intCvxComb) and with the direction taken from the
 * gradient or from the two-piece bundle (intFWDirection), has to find that
 * value, and the bound it reports has to hold on both sides of it. A last
 * case changes the costs of a sub-Block and solves again, which is what a
 * Solver attached to a Block has to survive.
 *
 * Two corners of the solver are in the instance on purpose: a coordinate
 * whose box is a single point, and a ColVariable of the father Objective that
 * the Objective of its sub-Block does not price, which FrankWolfeSolver adds
 * there with a zero cost to have somewhere to scatter the gradient, and takes
 * out again when it is unregistered.
 *
 * Another Solver registered to the father, which never computes, has to
 * receive nothing from a compute() of many iterations, which rewrites the
 * costs of all the sub-Block at each of them and sets them back at the end:
 * also when the father has a default channel of its own, which it has to
 * get back, and when the oracle of a sub-Block throws halfway, after which
 * the costs have to be the original ones and a new compute() has to find
 * the optimum.
 *
 * The test needs nothing but the core, so that the CI of FrankWolfeSolver
 * builds this module alone.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractBlock.h"
#include "BoxSolver.h"
#include "DQuadFunction.h"
#include "FakeSolver.h"
#include "FRealObjective.h"
#include "FrankWolfeSolver.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- DATA -------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

const int K = 3;           // number of sub-Block
const int N = 6;           // number of ColVariable of each sub-Block
const double tol = 1e-6;   // relative accuracy asked of the value

// the corners: the box of this coordinate is a single point, and the
// Objective of its sub-Block does not price this other one
const int fixed_coord = 2;
const int unpriced_coord = N + N - 1;

/// the data of one ColVariable: bounds, father coefficients, sub-Block cost
struct Coord {
 double l , u , a , b , c;
 };

std::vector< Coord > data;  // K * N entries, sub-Block by sub-Block

/// the optimal value of the whole problem, in closed form
double optimum( void )
{
 double v = 0;
 for( const auto & d : data ) {
  const double x = std::clamp( - ( d.b + d.c ) / ( 2 * d.a ) , d.l , d.u );
  v += d.a * x * x + ( d.b + d.c ) * x;
  }
 return( v );
 }

/*--------------------------------------------------------------------------*/
/// the father Block, its K sub-Block and the Solver of each of them

/// a BoxSolver whose compute() throws once it has been called limit times

class ThrowingBoxSolver : public BoxSolver {
 public:
 int compute( bool changedvars = true ) override {
  if( ++calls > limit )
   throw( std::runtime_error( "ThrowingBoxSolver::compute: limit reached" ) );
  return( BoxSolver::compute( changedvars ) );
  }
 int calls = 0;
 int limit = 1 << 30;
 };

AbstractBlock * build( std::vector< LinearFunction * > & costs ,
                       Solver * first = nullptr )
{
 auto father = new AbstractBlock();
 DQuadFunction::v_coeff_triple triples;
 triples.reserve( K * N );

 for( int j = 0 ; j < K ; ++j ) {
  auto sb = new AbstractBlock( father );

  auto x = new std::vector< ColVariable >( N );
  sb->add_static_variable( *x , "x" );

  auto box = new std::vector< BoxConstraint >( N );
  LinearFunction::v_coeff_pair pairs;
  for( int i = 0 ; i < N ; ++i ) {
   const auto & d = data[ j * N + i ];
   ( *box )[ i ].set_variable( & ( *x )[ i ] );
   ( *box )[ i ].set_lhs( d.l );
   ( *box )[ i ].set_rhs( d.u );
   if( j * N + i != unpriced_coord )
    pairs.emplace_back( & ( *x )[ i ] , d.c );
   triples.emplace_back( & ( *x )[ i ] , d.b , d.a );
   }
  sb->add_static_constraint( *box , "box" );

  auto lf = new LinearFunction( std::move( pairs ) );
  costs.push_back( lf );
  auto obj = new FRealObjective( sb , lf );
  obj->set_sense( Objective::eMin , eNoMod );
  sb->set_objective( obj );

  sb->register_Solver( ( first && ( j == 0 ) ) ? first
                                              : Solver::new_Solver( "BoxSolver" ) );
  father->add_nested_Block( sb );
  }

 auto obj = new FRealObjective( father ,
                                new DQuadFunction( std::move( triples ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 father->set_objective( obj );

 return( father );
 }

/*--------------------------------------------------------------------------*/
/// the parameters of one run of FrankWolfeSolver

struct Case {
 std::string name;
 int algorithm , cvx_comb , direction;
 };

void set_par( Solver * s , const std::string & name , int value )
{
 s->set_par( s->int_par_str2idx( name ) , value );
 }

void set_par( Solver * s , const std::string & name , double value )
{
 s->set_par( s->dbl_par_str2idx( name ) , value );
 }

/*--------------------------------------------------------------------------*/
/// solves and checks the value and the bound against the closed form

bool check( Solver * fw , const std::string & name )
{
 const double opt = optimum();
 const int status = fw->compute( false );
 const double value = fw->get_var_value();
 const double lb = fw->get_lb();
 const double ub = fw->get_ub();
 const double scale = std::max( 1.0 , std::abs( opt ) );

 // the value within the accuracy asked, and a bound that holds, i.e. the
 // lower one below the optimum and the upper one above it, up to that same
 // accuracy
 const bool ok = ( status >= Solver::kOK ) && ( status < Solver::kError ) &&
                 ( std::abs( value - opt ) <= tol * scale ) &&
                 ( lb <= opt + tol * scale ) && ( ub >= opt - tol * scale );

 std::cout << std::left << std::setw( 34 ) << name << std::right
           << std::scientific << std::setprecision( 9 )
           << " value " << value << " optimum " << opt
           << " lb " << lb << " status " << status
           << ( ok ? "  OK" : "  KO" ) << std::endl;
 return( ok );
 }

/*--------------------------------------------------------------------------*/
/// the linear coefficients of the sub-Block Objectives

std::vector< double > coefficients( const std::vector< LinearFunction * > & c )
{
 std::vector< double > v;
 for( auto lf : c )
  for( Block::Index i = 0 ; i < lf->get_num_active_var() ; ++i )
   v.push_back( lf->get_coefficient( i ) );
 return( v );
 }

/*--------------------------------------------------------------------------*/
/// what the other Solver of the father receives while FrankWolfeSolver runs

bool check_other_solver( void )
{
 bool ok = true;
 auto report = [ & ]( bool cond , const std::string & what ) {
  std::cout << std::left << std::setw( 60 ) << what
            << ( cond ? "  OK" : "  KO" ) << std::endl;
  ok &= cond;
  };

 std::vector< LinearFunction * > costs;
 auto thrower = new ThrowingBoxSolver();
 auto father = build( costs , thrower );
 auto fw = Solver::new_Solver( "FrankWolfeSolver" );
 father->register_Solver( fw );
 auto other = new FakeSolver();
 father->register_Solver( other );
 auto & mods = other->get_Modification_list();

 // vanilla with the gradient converges slowly, i.e., many iterations
 set_par( fw , "intLMOObj" , 2 );
 set_par( fw , "intAlgorithm" , 0 );
 set_par( fw , "intCvxComb" , 1 );
 set_par( fw , "intFWDirection" , 0 );
 set_par( fw , "intMaxIter" , 200000 );
 set_par( fw , "dblRelAcc" , 1e-10 );

 // the Variable FrankWolfeSolver adds to a sub-Block Objective when it is
 // registered are a real change, which the other Solver has already seen
 mods.clear();
 const auto c_before = coefficients( costs );

 ok &= check( fw , "with another Solver on the father" );
 const int iters = thrower->calls;
 report( iters > 100 , "  it took " + std::to_string( iters ) +
         " iterations" );
 report( mods.empty() , "  the other Solver has queued " +
         std::to_string( mods.size() ) + " Modification" );
 report( coefficients( costs ) == c_before , "  the costs are the original" );
 report( father->get_default_channel() == 0 ,
         "  the father has its default channel back" );

 // a default channel of the father set by someone else is given back, and
 // nothing is sent to it
 const auto ch = father->open_channel();
 father->set_default_channel( ch );
 mods.clear();
 ok &= check( fw , "with a default channel of the father" );
 report( father->get_default_channel() == ch ,
         "  the father has its default channel back" );
 father->close_channel( ch );
 bool empty = mods.size() <= 1;
 if( ! mods.empty() ) {
  auto gm = std::dynamic_pointer_cast< GroupModification >( mods.front() );
  empty &= gm && gm->sub_Modifications().empty();
  }
 report( empty , "  nothing was sent to it" );

 // a real change of the costs reaches the other Solver, and FrankWolfeSolver
 for( int i = 0 ; i < N ; ++i ) {
  data[ i ].c = - data[ i ].c;
  costs[ 0 ]->modify_coefficient( i , data[ i ].c );
  }
 report( ! mods.empty() , "  a change of the costs reaches the other Solver" );

 // the oracle of sub-Block 0 throws halfway: the exception comes out, the
 // costs are the original ones, nothing reached the other Solver, and the
 // next compute() finds the optimum
 mods.clear();
 const auto c_throw = coefficients( costs );
 thrower->calls = 0;
 thrower->limit = 5;
 bool threw = false;
 try {
  fw->compute( false );
  }
 catch( std::runtime_error & ) {
  threw = true;
  }
 report( threw , "  the exception of the oracle comes out of compute()" );
 report( mods.empty() , "  after it the other Solver has queued " +
         std::to_string( mods.size() ) + " Modification" );
 report( coefficients( costs ) == c_throw , "  after it the costs are the "
         "original" );
 report( father->get_default_channel() == 0 ,
         "  after it the father has its default channel back" );
 thrower->limit = 1 << 30;
 ok &= check( fw , "after the exception" );
 report( mods.empty() , "  the other Solver has queued " +
         std::to_string( mods.size() ) + " Modification" );

 // undo the change of the costs for the cases that may follow
 for( int i = 0 ; i < N ; ++i )
  data[ i ].c = - data[ i ].c;

 father->unregister_Solver( other );
 delete other;
 father->unregister_Solver( fw );
 delete fw;
 for( auto sb : father->get_nested_Blocks() ) {
  auto s = sb->get_registered_solvers().front();
  sb->unregister_Solver( s );
  delete s;
  }
 delete father;
 return( ok );
 }

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 // the data: some optimal coordinates in the interior of their box and some
 // on either of its faces, which the ranges below give with a fixed seed
 std::mt19937 rg( 7 );
 std::uniform_real_distribution< double > ul( -4 , 0 ) , uu( 0 , 4 ) ,
  ua( 0.5 , 2 ) , ub( -6 , 6 ) , uc( -3 , 3 );
 for( int k = 0 ; k < K * N ; ++k )
  data.push_back( { ul( rg ) , uu( rg ) , ua( rg ) , ub( rg ) , uc( rg ) } );
 data[ fixed_coord ].u = data[ fixed_coord ].l;
 data[ unpriced_coord ].c = 0;

 const std::vector< Case > cases = {
  { "vanilla, gradient" ,               0 , 1 , 0 } ,
  { "vanilla, bundle direction" ,       0 , 1 , 1 } ,
  { "Away-step, gradient" ,             1 , 1 , 0 } ,
  { "Away-step, bundle direction" ,     1 , 1 , 1 } ,
  { "BPCG, gradient" ,                  2 , 1 , 0 } ,
  { "BPCG, bundle direction" ,          2 , 1 , 1 } ,
  { "vanilla, objective at x" ,         0 , 0 , 0 } ,
  { "BPCG, objective at x" ,            2 , 0 , 0 } };

 bool all = true;
 for( const auto & c : cases ) {
  std::vector< LinearFunction * > costs;
  auto father = build( costs );
  auto fw = Solver::new_Solver( "FrankWolfeSolver" );
  father->register_Solver( fw );

  set_par( fw , "intLMOObj" , 2 );         // the sub-Block costs are in
  set_par( fw , "intAlgorithm" , c.algorithm );
  set_par( fw , "intCvxComb" , c.cvx_comb );
  set_par( fw , "intFWDirection" , c.direction );
  set_par( fw , "intMaxIter" , 200000 );
  set_par( fw , "dblRelAcc" , 1e-10 );

  all &= check( fw , c.name );

  // the last case also changes the costs of the first sub-Block and asks
  // for the new optimum, i.e. the Modification has to reach the Solver
  if( &c == &cases.back() ) {
   for( int i = 0 ; i < N ; ++i ) {
    data[ i ].c = - data[ i ].c;
    costs[ 0 ]->modify_coefficient( i , data[ i ].c );
    }
   all &= check( fw , c.name + ", costs changed" );
   }

  // unregistered, the solver has to give the sub-Block back its Objective
  // as it was, i.e. without the ColVariable it added there
  father->unregister_Solver( fw );
  delete fw;
  const auto n_after = costs[ 1 ]->get_num_active_var();
  if( n_after != Block::Index( N - 1 ) ) {
   std::cout << c.name << ": the Objective of sub-Block 1 has " << n_after
             << " ColVariable after the solver left, not " << N - 1
             << "  KO" << std::endl;
   all = false;
   }
  for( auto sb : father->get_nested_Blocks() ) {
   auto s = sb->get_registered_solvers().front();
   sb->unregister_Solver( s );
   delete s;
   }
  delete father;
  }

 all &= check_other_solver();

 std::cout << ( all ? "All tests passed!!" : "Some test FAILED!!" )
           << std::endl;
 return( all ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
