/*--------------------------------------------------------------------------*/
/*------------------ File IntegralityBarrierFunction.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntegralityBarrierFunction class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesca Demelas \n
 *         LIPN \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli, Francesca Demelas
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <unordered_map>

#include "Block.h"
#include "FRowConstraint.h"
#include "IntegralityBarrierFunction.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

using Index = Block::Index;
using Value = Function::FunctionValue;

// tightens the bounds lb[ i ] / ub[ i ] of the Variable with those of the
// OneVarConstraint of type C in the static group g of block
template< class C >
void bounds_of( Block * block , Index g ,
		const std::unordered_map< const Variable * , Index > & idx ,
		std::vector< Value > & lb , std::vector< Value > & ub )
{
 auto grp = block->get_static_constraint_v< C >( g );
 if( ! grp )
  return;
 for( const auto & c : *grp ) {
  auto it = idx.find( c.get_active_var( 0 ) );
  if( it == idx.end() )
   continue;
  lb[ it->second ] = std::max( lb[ it->second ] , Value( c.get_lhs() ) );
  ub[ it->second ] = std::min( ub[ it->second ] , Value( c.get_rhs() ) );
  }
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::build( Block * block )
{
 static const std::string _prfx = "IntegralityBarrierFunction::build: ";

 v_vars.clear();
 v_int.clear();
 v_rows.clear();

 // the Block and all those inside it, a father before its sub-Block
 std::vector< Block * > tree = { block };
 for( std::size_t t = 0 ; t < tree.size() ; ++t )
  for( auto sub : tree[ t ]->get_nested_Blocks() )
   tree.push_back( sub );

 // the Variable - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 std::unordered_map< const Variable * , Index > idx;
 auto add_var = [ & ]( ColVariable & x ) {
  idx.emplace( & x , v_vars.size() );
  v_vars.push_back( & x );
  v_int.push_back( x.is_integer() );
  };
 for( auto b : tree ) {
  for( Index g = 0 ; g < b->get_number_static_variables() ; ++g )
   if( auto grp = b->get_static_variable_v< ColVariable >( g ) )
    for( auto & x : *grp )
     add_var( x );
  for( Index g = 0 ; g < b->get_number_dynamic_variables() ; ++g )
   if( auto grp = b->get_dynamic_variable_v< ColVariable >( g ) )
    for( auto & lst : *grp )
     for( auto & x : lst )
      add_var( x );
  }
 const Index n = v_vars.size();

 // their bounds - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 std::vector< Value > lb( n ) , ub( n );
 for( Index i = 0 ; i < n ; ++i ) {
  lb[ i ] = v_vars[ i ]->get_lb();
  ub[ i ] = v_vars[ i ]->get_ub();
  }
 for( auto b : tree )
  for( Index g = 0 ; g < b->get_number_static_constraints() ; ++g ) {
   bounds_of< BoxConstraint >( b , g , idx , lb , ub );
   bounds_of< LBConstraint >( b , g , idx , lb , ub );
   bounds_of< UBConstraint >( b , g , idx , lb , ub );
   bounds_of< LB0Constraint >( b , g , idx , lb , ub );
   bounds_of< UB0Constraint >( b , g , idx , lb , ub );
   bounds_of< NNConstraint >( b , g , idx , lb , ub );
   bounds_of< NPConstraint >( b , g , idx , lb , ub );
   bounds_of< ZOConstraint >( b , g , idx , lb , ub );
   }

 // the rows a x >= b, the equalities left out - - - - - - - - - - - - - - -
 auto add_row = [ & ]( const LinearFunction::v_coeff_pair & cp , Value sign ,
		       Value rhs ) {
  Row r;
  r.rhs = sign * rhs;
  Value maxact = 0;  // the largest a x over the bounds
  for( const auto & [ var , a ] : cp ) {
   if( a == 0 )
    continue;
   auto it = idx.find( var );
   if( it == idx.end() )
    throw( std::invalid_argument( _prfx + "a row has a Variable that is "
				  "not a ColVariable of the Block" ) );
   r.idx.push_back( it->second );
   r.coef.push_back( sign * a );
   maxact += sign * a * ( sign * a > 0 ? ub[ it->second ] : lb[ it->second ] );
   }
  // the largest slack the row can have, if finite and positive
  const Value range = maxact - r.rhs;
  r.norm = ( std::isfinite( range ) && ( range > 0 ) ) ? range : Value( 1 );
  v_rows.push_back( std::move( r ) );
  };

 auto add_rows_of = [ & ]( FRowConstraint & c ) {
   auto lf = dynamic_cast< LinearFunction * >( c.get_function() );
   if( ! lf )
    throw( std::invalid_argument( _prfx + "a FRowConstraint whose Function "
				  "is not a LinearFunction" ) );
   const Value lhs = c.get_lhs();
   const Value rhs = c.get_rhs();
   if( c.is_relaxed() || ( lhs == rhs ) )  // a relaxed row or an equality,
    return;                             // whose slack is always 0
   const Value cst = lf->get_constant_term();
   if( lhs > - Inf< Value >() )
    add_row( lf->get_v_var() , 1 , lhs - cst );
   if( rhs < Inf< Value >() )
    add_row( lf->get_v_var() , -1 , rhs - cst );
  };
 for( auto b : tree ) {
  for( Index g = 0 ; g < b->get_number_static_constraints() ; ++g )
   if( auto grp = b->get_static_constraint_v< FRowConstraint >( g ) )
    for( auto & c : *grp )
     add_rows_of( c );
  for( Index g = 0 ; g < b->get_number_dynamic_constraints() ; ++g )
   if( auto grp = b->get_dynamic_constraint_v< FRowConstraint >( g ) )
    for( auto & lst : *grp )
     for( auto & c : lst )
      add_rows_of( c );
  }

 v_rows_of.assign( n , {} );
 for( Index k = 0 ; k < v_rows.size() ; ++k )
  for( auto i : v_rows[ k ].idx )
   v_rows_of[ i ].push_back( k );

 v_y.assign( v_rows.size() , 1 );
 v_grad.assign( n , 0 );
 v_x.assign( n , 0 );
 f_value = 0;
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::set_y( std::vector< FunctionValue > && y ,
					ModParam issueMod )
{
 if( y.size() != v_rows.size() )
  throw( std::invalid_argument( "IntegralityBarrierFunction::set_y: " +
				std::to_string( y.size() ) + " exponents for " +
				std::to_string( v_rows.size() ) + " rows" ) );
 for( auto yk : y )
  if( ! ( yk >= 0 ) )
   throw( std::invalid_argument( "IntegralityBarrierFunction::set_y: a "
				 "negative exponent" ) );
 v_y = std::move( y );
 everything_changed( issueMod );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::set_eps( FunctionValue eps ,
					  ModParam issueMod )
{
 if( ! ( eps > 0 ) )
  throw( std::invalid_argument( "IntegralityBarrierFunction::set_eps: "
				"epsilon must be positive" ) );
 f_eps = eps;
 everything_changed( issueMod );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::set_psi( int psi , ModParam issueMod )
{
 if( ( psi != ePsiPlain ) && ( psi != ePsiSqrt ) )
  throw( std::invalid_argument( "IntegralityBarrierFunction::set_psi: "
				"unknown form " + std::to_string( psi ) ) );
 f_psi = psi;
 everything_changed( issueMod );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::everything_changed( ModParam issueMod )
{
 if( ( ! f_Observer ) || ( ! f_Observer->issue_mod( issueMod ) ) )
  return;
 f_Observer->add_Modification( std::make_shared< FunctionMod >( this ) ,
			       Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/
/*---------- METHODS DESCRIBING THE BEHAVIOR OF THE C05Function ------------*/
/*--------------------------------------------------------------------------*/

Function::FunctionValue IntegralityBarrierFunction::psi( FunctionValue t )
 const
{
 FunctionValue B = std::cos( 2 * M_PI * ( t - 0.5 ) ) + 1 + f_eps;
 if( f_psi == ePsiSqrt )
  B = std::sqrt( B );
 return( std::exp( - 1 / B ) );
 }

/*--------------------------------------------------------------------------*/

Function::FunctionValue IntegralityBarrierFunction::dpsi( FunctionValue t )
 const
{
 // psi = exp( - 1 / B ), B = cos( 2 pi ( t - 1/2 ) ) + 1 + eps, hence
 // psi' = psi * B' / B^2 with B' = - 2 pi sin( 2 pi ( t - 1/2 ) ); with
 // ePsiSqrt B is the square root of that, and B' is divided by 2 B
 FunctionValue B = std::cos( 2 * M_PI * ( t - 0.5 ) ) + 1 + f_eps;
 FunctionValue dB = - 2 * M_PI * std::sin( 2 * M_PI * ( t - 0.5 ) );
 if( f_psi == ePsiSqrt ) {
  B = std::sqrt( B );
  dB /= 2 * B;
  }
 return( std::exp( - 1 / B ) * dB / ( B * B ) );
 }

/*--------------------------------------------------------------------------*/

Function::FunctionValue IntegralityBarrierFunction::slack( Index k ) const
{
 const auto & r = v_rows[ k ];
 FunctionValue ax = 0;
 for( Index p = 0 ; p < r.idx.size() ; ++p )
  ax += r.coef[ p ] * v_vars[ r.idx[ p ] ]->get_value();
 return( ( ax - r.rhs ) / r.norm );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::barrier(
				       std::vector< FunctionValue > & sigma ,
				       std::vector< bool > & unclamped ,
				       std::vector< FunctionValue > & P ) const
{
 const Index m = v_rows.size();
 sigma.resize( m );
 unclamped.resize( m );
 for( Index k = 0 ; k < m ; ++k ) {
  const FunctionValue t = std::max( slack( k ) , FunctionValue( 0 ) ) + f_eps;
  unclamped[ k ] = ( t < 1 );
  sigma[ k ] = std::min( t , FunctionValue( 1 ) );
  }

 P.assign( v_vars.size() , 1 );
 for( Index i = 0 ; i < v_vars.size() ; ++i )
  if( v_int[ i ] )
   for( auto k : v_rows_of[ i ] )
    if( v_y[ k ] != 0 )
     P[ i ] *= std::pow( sigma[ k ] , v_y[ k ] );
 }

/*--------------------------------------------------------------------------*/

int IntegralityBarrierFunction::compute( bool changedvars )
{
 const Index n = v_vars.size();
 std::vector< FunctionValue > sigma , P;
 std::vector< bool > unclamped;
 barrier( sigma , unclamped , P );

 /* f = sum_i psi_i / D_i with D_i = P_i + eps over the integer i, hence
  *
  *  df / dx_j = psi'_j / D_j - sum_i psi_i / D_i^2 dP_i / dx_j
  *
  * and dP_i / dx_j = P_i sum_{k in C( i )} ( y_k / sigma_k ) dsigma_k / dx_j,
  * where dsigma_k / dx_j = a_kj / n_k if sigma_k is not cut at 1, and 0
  * otherwise; so the second term is A^T r, with
  *
  *  r_k = ( y_k / ( sigma_k n_k ) ) sum_{i : k in C( i )} q_i ,
  *  q_i = - psi_i P_i / D_i^2 . */

 f_value = 0;
 v_grad.assign( n , 0 );
 std::vector< FunctionValue > q( n , 0 );
 for( Index i = 0 ; i < n ; ++i ) {
  v_x[ i ] = v_vars[ i ]->get_value();
  if( ! v_int[ i ] )
   continue;
  const FunctionValue ps = psi( v_x[ i ] );
  const FunctionValue D = P[ i ] + f_eps;
  f_value += ps / D;
  v_grad[ i ] += dpsi( v_x[ i ] ) / D;
  q[ i ] = - ps * P[ i ] / ( D * D );
  }

 for( Index k = 0 ; k < v_rows.size() ; ++k ) {
  if( ( v_y[ k ] == 0 ) || ( ! unclamped[ k ] ) )
   continue;
  const auto & r = v_rows[ k ];
  FunctionValue sq = 0;
  for( auto i : r.idx )
   if( v_int[ i ] )
    sq += q[ i ];
  if( sq == 0 )
   continue;
  const FunctionValue rk = v_y[ k ] * sq / ( sigma[ k ] * r.norm );
  for( Index p = 0 ; p < r.idx.size() ; ++p )
   v_grad[ r.idx[ p ] ] += r.coef[ p ] * rk;
  }

 return( kOK );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::get_y_gradient(
				    std::vector< FunctionValue > & gy ) const
{
 std::vector< FunctionValue > sigma , P;
 std::vector< bool > unclamped;
 barrier( sigma , unclamped , P );

 gy.assign( v_rows.size() , 0 );
 for( Index i = 0 ; i < v_vars.size() ; ++i ) {
  if( ! v_int[ i ] )
   continue;
  const FunctionValue D = P[ i ] + f_eps;
  const FunctionValue w = - psi( v_vars[ i ]->get_value() ) * P[ i ] /
			  ( D * D );
  for( auto k : v_rows_of[ i ] )
   gy[ k ] += w * std::log( sigma[ k ] );
  }
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::get_linearization_coefficients(
			       FunctionValue * g , Range range , Index name )
{
 if( name != Inf< Index >() )
  throw( std::logic_error( "IntegralityBarrierFunction::"
			   "get_linearization_coefficients: only the "
			   "diagonal linearization is there" ) );
 range.second = std::min( range.second , Index( v_grad.size() ) );
 for( Index i = range.first ; i < range.second ; ++i )
  g[ i - range.first ] = v_grad[ i ];
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::get_linearization_coefficients(
			     FunctionValue * g , c_Subset & subset ,
			     bool ordered , Index name )
{
 if( name != Inf< Index >() )
  throw( std::logic_error( "IntegralityBarrierFunction::"
			   "get_linearization_coefficients: only the "
			   "diagonal linearization is there" ) );
 for( auto i : subset )
  *(g++) = v_grad[ i ];
 }

/*--------------------------------------------------------------------------*/

Function::FunctionValue IntegralityBarrierFunction::get_linearization_constant(
								  Index name )
{
 if( name != Inf< Index >() )
  throw( std::logic_error( "IntegralityBarrierFunction::"
			   "get_linearization_constant: only the diagonal "
			   "linearization is there" ) );
 FunctionValue c = f_value;
 for( Index i = 0 ; i < v_grad.size() ; ++i )
  c -= v_grad[ i ] * v_x[ i ];
 return( c );
 }

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING "ACTIVE" Variable IN THE C05Function --------*/
/*--------------------------------------------------------------------------*/

Block::Index IntegralityBarrierFunction::is_active( const Variable * var )
 const
{
 for( Index i = 0 ; i < v_vars.size() ; ++i )
  if( v_vars[ i ] == var )
   return( i );
 return( Inf< Index >() );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierFunction::print( std::ostream & output )
{
 output << "IntegralityBarrierFunction with " << v_vars.size()
	<< " Variable and " << v_rows.size() << " rows" << std::endl;
 }

/*--------------------------------------------------------------------------*/
/*---------------- End File IntegralityBarrierFunction.cpp -----------------*/
/*--------------------------------------------------------------------------*/
