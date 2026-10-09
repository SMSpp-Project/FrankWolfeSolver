/*--------------------------------------------------------------------------*/
/*------------------- File IntegralityBarrierFunction.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the IntegralityBarrierFunction class, a C05Function whose
 * local minima over the polyhedron of a (mixed-)integer linear program are,
 * hopefully, only its integer points: an integrality penalty, zero exactly
 * on the integer values, divided by a barrier-like product of the slacks of
 * the constraints, which is small on the boundary of the polyhedron, so that
 * the non-integer points of the boundary are penalised far more than those
 * in the interior, while the integer ones, the penalty being zero there,
 * stay global minima even if they lie on the boundary.
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
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli, Francesca Demelas
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IntegralityBarrierFunction
 #define __IntegralityBarrierFunction
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "C05Function.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------- CLASS IntegralityBarrierFunction --------------------*/
/*--------------------------------------------------------------------------*/
/// an integrality penalty divided by a product of slacks
/** The IntegralityBarrierFunction is
 *
 *  \f[ f( x ) = \sum_{i \in I} \frac{ \psi( x_i ) }
 *                                   { P_i( x ) + \varepsilon } \;\; , \;\;
 *      P_i( x ) = \prod_{k \in C( i )} \sigma_k( x )^{y_k} \f]
 *
 * where I is the set of the integer variables, C( i ) that of the rows in
 * which x_i appears, and
 *
 *  \f[ \psi( t ) = e^{ - 1 / ( \cos( 2 \pi ( t - 1/2 ) ) + 1 + \varepsilon ) }
 *  \f]
 *
 * or, with ePsiSqrt [see set_psi()],
 *
 *  \f[ \psi( t ) = e^{ - 1 / \sqrt{ \cos( 2 \pi ( t - 1/2 ) ) + 1 +
 *                                     \varepsilon } } \f]
 *
 * is a smooth integrality penalty, periodic of period 1, smallest on the
 * integer values and largest on the half-integer ones; the square root
 * makes it go to 0 faster near the integer values. The rows are those of
 * the polyhedron written as \f$ a_k x \geq b_k \f$, each one with its slack
 * \f$ s_k = ( a_k x - b_k ) / n_k \f$, normalised by the largest value
 * \f$ n_k \f$ the slack can take over the bounds of the variables (1 if that
 * is not finite and positive), and
 *
 *  \f[ \sigma_k( x ) = \min\{ \max\{ s_k , 0 \} + \varepsilon \,,\, 1 \} \f]
 *
 * so that the barrier is the product of the normalised slacks, kept away
 * from 0 by \f$ \varepsilon \f$ and from growing past 1. The exponents
 * \f$ y_k \geq 0 \f$ weigh the rows: with \f$ y_k = 0 \f$ the slack of row k
 * does not count at all, with \f$ y_k = 1 \f$ it counts as it is [see
 * set_y()]; they all start at 1.
 *
 * The rows are made out of a Block [see build()]: those of its
 * FRowConstraint, each one a LinearFunction, with a finite lhs giving
 * \f$ a x \geq lhs \f$ and a finite rhs giving \f$ - a x \geq - rhs \f$,
 * while the equalities (lhs == rhs) and the bounds of the variables are left
 * out, since their slack is always 0 or is taken care of by \f$ \psi \f$.
 *
 * The function is differentiable, and compute() computes its value and its
 * gradient, which is the (diagonal) linearization it gives; the gradient
 * with respect to the exponents y is given by get_y_gradient(). */

class IntegralityBarrierFunction : public C05Function
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 using v_col_var = std::vector< ColVariable * >;  ///< the active Variable

 /// the forms of the integrality penalty psi [see the class]
 enum psi_type {
  ePsiPlain = 0 ,  ///< e^{ - 1 / ( cos + 1 + eps ) }
  ePsiSqrt = 1     ///< e^{ - 1 / sqrt( cos + 1 + eps ) }
  };

 /// a row a x >= b: the indices of its nonzeros, their values, and b
 struct Row {
  std::vector< Index > idx;            ///< the active Variable in the row
  std::vector< FunctionValue > coef;   ///< their coefficients
  FunctionValue rhs;                   ///< the b of a x >= b
  FunctionValue norm;                  ///< the normalisation of the slack
  };

/*--------------------------------------------------------------------------*/
 /// virtualized concrete iterator
 class v_iterator : public ThinVarDepInterface::v_iterator
 {
  public:
  explicit v_iterator( v_col_var::iterator itr ) : itr_( itr ) {}
  v_iterator * clone( void ) override { return( new v_iterator( itr_ ) ); }
  void operator++( void ) override final { ++itr_; }
  reference operator*( void ) const override final { return( **itr_ ); }
  pointer operator->( void ) const override final { return( *itr_ ); }
  bool operator==( const ThinVarDepInterface::v_iterator & rhs )
   const override final {
   auto tmp = dynamic_cast< const v_iterator * >( & rhs );
   return( tmp ? itr_ == tmp->itr_ : false );
   }
  bool operator!=( const ThinVarDepInterface::v_iterator & rhs )
   const override final {
   auto tmp = dynamic_cast< const v_iterator * >( & rhs );
   return( tmp ? itr_ != tmp->itr_ : true );
   }
  private:
  v_col_var::iterator itr_;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// virtualized concrete const_iterator
 class v_const_iterator : public ThinVarDepInterface::v_const_iterator
 {
  public:
  explicit v_const_iterator( v_col_var::const_iterator itr ) : itr_( itr ) {}
  v_const_iterator * clone( void ) override {
   return( new v_const_iterator( itr_ ) );
   }
  void operator++( void ) override final { ++itr_; }
  reference operator*( void ) const override final { return( **itr_ ); }
  pointer operator->( void ) const override final { return( *itr_ ); }
  bool operator==( const ThinVarDepInterface::v_const_iterator & rhs )
   const override final {
   auto tmp = dynamic_cast< const v_const_iterator * >( & rhs );
   return( tmp ? itr_ == tmp->itr_ : false );
   }
  bool operator!=( const ThinVarDepInterface::v_const_iterator & rhs )
   const override final {
   auto tmp = dynamic_cast< const v_const_iterator * >( & rhs );
   return( tmp ? itr_ != tmp->itr_ : true );
   }
  private:
  v_col_var::const_iterator itr_;
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: an empty function, to be built() out of a Block
 /** Constructor: the function has no Variable and no row until build() is
  * called; \p eps is the \f$ \varepsilon \f$ of the formulae [see the
  * class]. */

 explicit IntegralityBarrierFunction( FunctionValue eps = 1e-6 )
  : C05Function() , f_eps( eps ) , f_value( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor: nothing to do

 ~IntegralityBarrierFunction() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /// makes the function out of the Variable and the rows of a Block
 /** Makes the function out of the Block \p block: its active Variable are
  * the ColVariable in the static and then the dynamic groups of \p block
  * and of all the Block inside it, a father before its sub-Block, in their
  * order, and its rows those of the FRowConstraint, but the relaxed ones, in
  * the same groups [see the class], whose LinearFunction must have no
  * Variable but those. The bounds of the Variable, needed by the
  * normalisation of the slacks, are those of the ColVariable together with
  * those of the OneVarConstraint in the static groups. The exponents y are
  * all set to 1. The function issues no Modification: it is meant to be
  * built before it is used. */

 void build( Block * block );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the exponents y of the rows
 /** Sets the exponents y of the rows, one per row [see get_rows()], each
  * one non-negative; a FunctionMod saying that everything has changed is
  * issued as \p issueMod says. */

 void set_y( std::vector< FunctionValue > && y ,
	     ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the epsilon of the formulae

 void set_eps( FunctionValue eps , ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the form of the integrality penalty, a psi_type value

 void set_psi( int psi , ModParam issueMod = eModBlck );

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
 /// removing Variable is not supported

 void remove_variable( Index i , ModParam issueMod = eModBlck ) override {
  throw( std::logic_error( "IntegralityBarrierFunction::remove_variable: "
			   "not supported" ) );
  }

/*--------------------------------------------------------------------------*/
/*---------- METHODS DESCRIBING THE BEHAVIOR OF THE C05Function ------------*/
/*--------------------------------------------------------------------------*/
 /// computes the value and the gradient at the values of the Variable

 int compute( bool changedvars = true ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the value computed by the last compute()

 FunctionValue get_value( void ) override { return( f_value ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the function is not convex

 bool is_convex( void ) override { return( false ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the function is continuously differentiable

 bool is_continuously_differentiable( void ) const override {
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the coefficients of the gradient of the last compute(), in a Range

 void get_linearization_coefficients( FunctionValue * g ,
				      Range range = INFRange ,
				      Index name = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the coefficients of the gradient of the last compute(), in a Subset

 void get_linearization_coefficients( FunctionValue * g , c_Subset & subset ,
				      bool ordered = false ,
				      Index name = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the constant of the linearization of the last compute()
 /** The constant \f$ f( \bar{x} ) - g \bar{x} \f$ of the linearization at
  * the point \f$ \bar{x} \f$ of the last compute(). */

 FunctionValue get_linearization_constant( Index name = Inf< Index >() )
  override;

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE FUNCTION --------------*/
/*--------------------------------------------------------------------------*/
 /// the rows a x >= b of the barrier

 [[nodiscard]] const std::vector< Row > & get_rows( void ) const {
  return( v_rows );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the exponents y of the rows

 [[nodiscard]] const std::vector< FunctionValue > & get_y( void ) const {
  return( v_y );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the normalised slack of row k at the values of the Variable

 [[nodiscard]] FunctionValue slack( Index k ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the gradient with respect to y at the values of the Variable
 /** Writes in \p gy the partial derivatives of the function with respect to
  * the exponents y, at the current values of the Variable:
  *
  *  \f[ \partial f / \partial y_k = - \sum_{i \in I : k \in C( i )}
  *      \frac{ \psi( x_i ) P_i( x ) }{ ( P_i( x ) + \varepsilon )^2 }
  *      \ln \sigma_k( x ) \f] */

 void get_y_gradient( std::vector< FunctionValue > & gy ) const;

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING "ACTIVE" Variable IN THE C05Function --------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] Index get_num_active_var( void ) const override {
  return( v_vars.size() );
  }

 [[nodiscard]] Index is_active( const Variable * var ) const override;

 [[nodiscard]] Variable * get_active_var( Index i ) const override {
  return( v_vars[ i ] );
  }

 v_iterator * v_begin( void ) override {
  return( new v_iterator( v_vars.begin() ) );
  }

 [[nodiscard]] v_const_iterator * v_begin( void ) const override {
  return( new v_const_iterator( v_vars.cbegin() ) );
  }

 v_iterator * v_end( void ) override {
  return( new v_iterator( v_vars.end() ) );
  }

 [[nodiscard]] v_const_iterator * v_end( void ) const override {
  return( new v_const_iterator( v_vars.cend() ) );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
 /// prints the function: the number of Variable and of rows

 void print( std::ostream & output ) override;

/*--------------------------------------------------------------------------*/
 /// psi( t ) and its derivative

 [[nodiscard]] FunctionValue psi( FunctionValue t ) const;

 [[nodiscard]] FunctionValue dpsi( FunctionValue t ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the sigma of every row and the P of every integer Variable
 /** Computes, at the current values of the Variable, sigma[ k ] for every
  * row, unclamped[ k ] true if it is not cut at 1 (hence it depends on x),
  * and P[ i ] for every Variable (1 for the continuous ones). */

 void barrier( std::vector< FunctionValue > & sigma ,
	       std::vector< bool > & unclamped ,
	       std::vector< FunctionValue > & P ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// issues a FunctionMod saying that everything has changed

 void everything_changed( ModParam issueMod );

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 FunctionValue f_eps;                  ///< the epsilon of the formulae

 int f_psi = ePsiPlain;                ///< the form of psi

 v_col_var v_vars;                     ///< the active Variable

 std::vector< bool > v_int;            ///< which of them are integer

 std::vector< Row > v_rows;            ///< the rows a x >= b

 std::vector< std::vector< Index > > v_rows_of;  ///< C( i ) for each i

 std::vector< FunctionValue > v_y;     ///< the exponents of the rows

 FunctionValue f_value;                ///< the value of the last compute()

 std::vector< FunctionValue > v_grad;  ///< the gradient of the last compute()

 std::vector< FunctionValue > v_x;     ///< the point of the last compute()

/*--------------------------------------------------------------------------*/

 };  // end( class( IntegralityBarrierFunction ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* IntegralityBarrierFunction.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File IntegralityBarrierFunction.h ---------------*/
/*--------------------------------------------------------------------------*/
