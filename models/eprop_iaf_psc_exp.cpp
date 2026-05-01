/*
 *  eprop_iaf_psc_exp.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// nest models
#include "eprop_iaf_psc_exp.h"

// C++
#include <limits>

// libnestutil
#include "dict_util.h"
#include "iaf_propagator.h"
#include "numerics.h"

// nestkernel
#include "eprop_archiving_node_recurrent_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_impl.h"
#include "ring_buffer_impl.h"
#include "universal_data_logger_impl.h"

namespace nest
{

void
register_eprop_iaf_psc_exp( const std::string& name )
{
  register_node_model< eprop_iaf_psc_exp >( name );
}

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

RecordablesMap< eprop_iaf_psc_exp > eprop_iaf_psc_exp::recordablesMap_;

template <>
void
RecordablesMap< eprop_iaf_psc_exp >::create()
{
  // use standard names wherever you can for consistency!
  insert_( names::eprop_history_duration, &eprop_iaf_psc_exp::get_eprop_history_duration );
  insert_( names::V_m, &eprop_iaf_psc_exp::get_V_m_ );
  insert_( names::I_syn_ex, &eprop_iaf_psc_exp::get_I_syn_ex_ );
  insert_( names::I_syn_in, &eprop_iaf_psc_exp::get_I_syn_in_ );
  insert_( names::learning_signal, &eprop_iaf_psc_exp::get_learning_signal_ );
  insert_( names::surrogate_gradient, &eprop_iaf_psc_exp::get_surrogate_gradient_ );
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

eprop_iaf_psc_exp::Parameters_::Parameters_()
  : Tau_( 10.0 )              // in ms
  , C_( 250.0 )               // in pF
  , t_ref_( 2.0 )             // in ms
  , E_L_( -70.0 )             // in mV
  , I_e_( 0.0 )               // in pA
  , Theta_( -55.0 - E_L_ )    // relative E_L_
  , V_reset_( -70.0 - E_L_ )  // in mV
  , tau_ex_( 2.0 )            // in ms
  , tau_in_( 2.0 )            // in ms
  , c_reg_( 0.0 )
  , f_target_( 0.01 )
  , beta_( 1.0 )
  , gamma_( 0.3 )
  , surrogate_gradient_function_( "piecewise_linear" )
  , kappa_( 0.97 )
  , kappa_reg_( 0.97 )
  , eprop_isi_trace_cutoff_( 1000.0 )
{
}

eprop_iaf_psc_exp::State_::State_()
  : i_0_( 0.0 )
  , i_1_( 0.0 )
  , i_syn_ex_( 0.0 )
  , i_syn_in_( 0.0 )
  , V_m_( 0.0 )
  , r_ref_( 0 )
  , learning_signal_( 0.0 )
  , surrogate_gradient_( 0.0 )
{
}

eprop_iaf_psc_exp::Buffers_::Buffers_( eprop_iaf_psc_exp& n )
  : logger_( n )
{
}

eprop_iaf_psc_exp::Buffers_::Buffers_( const Buffers_&, eprop_iaf_psc_exp& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
eprop_iaf_psc_exp::Parameters_::get( Dictionary& d ) const
{
  d[ names::E_L ] = E_L_;  // resting potential
  d[ names::I_e ] = I_e_;
  d[ names::V_th ] = Theta_ + E_L_;  // threshold value
  d[ names::V_reset ] = V_reset_ + E_L_;
  d[ names::C_m ] = C_;
  d[ names::tau_m ] = Tau_;
  d[ names::tau_syn_ex ] = tau_ex_;
  d[ names::tau_syn_in ] = tau_in_;
  d[ names::t_ref ] = t_ref_;
  d[ names::c_reg ] = c_reg_;
  d[ names::f_target ] = f_target_;
  d[ names::beta ] = beta_;
  d[ names::gamma ] = gamma_;
  d[ names::surrogate_gradient_function ] = surrogate_gradient_function_;
  d[ names::kappa ] = kappa_;
  d[ names::kappa_reg ] = kappa_reg_;
  d[ names::eprop_isi_trace_cutoff ] = eprop_isi_trace_cutoff_;
}

double
eprop_iaf_psc_exp::Parameters_::set( const Dictionary& d, Node* node )
{
  // if E_L_ is changed, we need to adjust all variables defined relative to E_L_
  const double ELold = E_L_;
  update_value_param( d, names::E_L, E_L_, node );
  const double delta_EL = E_L_ - ELold;

  if ( update_value_param( d, names::V_reset, V_reset_, node ) )
  {
    V_reset_ -= E_L_;
  }
  else
  {
    V_reset_ -= delta_EL;
  }

  if ( update_value_param( d, names::V_th, Theta_, node ) )
  {
    Theta_ -= E_L_;
  }
  else
  {
    Theta_ -= delta_EL;
  }

  update_value_param( d, names::I_e, I_e_, node );
  update_value_param( d, names::C_m, C_, node );
  update_value_param( d, names::tau_m, Tau_, node );
  update_value_param( d, names::tau_syn_ex, tau_ex_, node );
  update_value_param( d, names::tau_syn_in, tau_in_, node );
  update_value_param( d, names::t_ref, t_ref_, node );

  update_value_param( d, names::c_reg, c_reg_, node );

  if ( update_value_param( d, names::f_target, f_target_, node ) )
  {
    f_target_ /= 1000.0;  // convert from spikes/s to spikes/ms
  }

  update_value_param( d, names::beta, beta_, node );
  update_value_param( d, names::gamma, gamma_, node );

  if ( update_value_param( d, names::surrogate_gradient_function, surrogate_gradient_function_, node ) )
  {
    eprop_iaf_psc_exp* nrn = dynamic_cast< eprop_iaf_psc_exp* >( node );
    assert( nrn );
    auto compute_surrogate_gradient = nrn->find_surrogate_gradient( surrogate_gradient_function_ );
    nrn->compute_surrogate_gradient_ = compute_surrogate_gradient;
  }

  update_value_param( d, names::kappa, kappa_, node );
  update_value_param( d, names::kappa_reg, kappa_reg_, node );
  update_value_param( d, names::eprop_isi_trace_cutoff, eprop_isi_trace_cutoff_, node );

  if ( V_reset_ >= Theta_ )
  {
    throw BadProperty( "Reset potential must be smaller than threshold." );
  }

  if ( C_ <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }

  if ( Tau_ <= 0 or tau_ex_ <= 0 or tau_in_ <= 0 )
  {
    throw BadProperty( "Membrane and synapse time constants must be strictly positive." );
  }

  if ( t_ref_ < 0 )
  {
    throw BadProperty( "Refractory time must not be negative." );
  }

  if ( c_reg_ < 0 )
  {
    throw BadProperty( "Firing rate regularization coefficient c_reg ≥ 0 required." );
  }

  if ( f_target_ < 0 )
  {
    throw BadProperty( "Firing rate regularization target rate f_target ≥ 0 required." );
  }

  if ( kappa_ < 0.0 or kappa_ > 1.0 )
  {
    throw BadProperty( "Eligibility trace low-pass filter kappa from range [0, 1] required." );
  }

  if ( kappa_reg_ < 0.0 or kappa_reg_ > 1.0 )
  {
    throw BadProperty( "Firing rate low-pass filter for regularization kappa_reg from range [0, 1] required." );
  }

  if ( eprop_isi_trace_cutoff_ < 0.0 )
  {
    throw BadProperty( "Cutoff of integration of eprop trace between spikes eprop_isi_trace_cutoff ≥ 0 required." );
  }

  return delta_EL;
}

void
eprop_iaf_psc_exp::State_::get( Dictionary& d, const Parameters_& p ) const
{
  d[ names::V_m ] = V_m_ + p.E_L_;  // Membrane potential
  d[ names::surrogate_gradient ] = surrogate_gradient_;
  d[ names::learning_signal ] = learning_signal_;
}

void
eprop_iaf_psc_exp::State_::set( const Dictionary& d, const Parameters_& p, double delta_EL, Node* node )
{
  if ( update_value_param( d, names::V_m, V_m_, node ) )
  {
    V_m_ -= p.E_L_;
  }
  else
  {
    V_m_ -= delta_EL;
  }
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

eprop_iaf_psc_exp::eprop_iaf_psc_exp()
  : EpropArchivingNodeRecurrent()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

eprop_iaf_psc_exp::eprop_iaf_psc_exp( const eprop_iaf_psc_exp& n )
  : EpropArchivingNodeRecurrent( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
eprop_iaf_psc_exp::init_buffers_()
{
  B_.input_buffer_.clear();  // includes resize
  B_.logger_.reset();
}

void
eprop_iaf_psc_exp::pre_run_hook()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  V_.eprop_isi_trace_cutoff_steps_ = Time( Time::ms( P_.eprop_isi_trace_cutoff_ ) ).get_steps();

  const double h = Time::get_resolution().get_ms();

  // these P are independent
  V_.P11ex_ = std::exp( -h / P_.tau_ex_ );
  V_.P11in_ = std::exp( -h / P_.tau_in_ );

  V_.P22_ = std::exp( -h / P_.Tau_ );

  // these are determined according to a numeric stability criterion
  V_.P21ex_ = IAFPropagatorExp( P_.tau_ex_, P_.Tau_, P_.C_ ).evaluate( h );
  V_.P21in_ = IAFPropagatorExp( P_.tau_in_, P_.Tau_, P_.C_ ).evaluate( h );

  V_.P20_ = P_.Tau_ / P_.C_ * ( 1.0 - V_.P22_ );

  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );
}

/* ----------------------------------------------------------------
 * Update function
 * ---------------------------------------------------------------- */

void
eprop_iaf_psc_exp::update( const Time& origin, const long from, const long to )
{
  for ( long lag = from; lag < to; ++lag )
  {
    const long t = origin.get_steps() + lag;

    if ( S_.r_ref_ == 0 )  // neuron not refractory, so evolve V
    {
      S_.V_m_ =
        S_.V_m_ * V_.P22_ + S_.i_syn_ex_ * V_.P21ex_ + S_.i_syn_in_ * V_.P21in_ + ( P_.I_e_ + S_.i_0_ ) * V_.P20_;
    }
    else
    {
      // neuron is absolute refractory
      --S_.r_ref_;
    }

    // exponential decaying PSCs
    S_.i_syn_ex_ *= V_.P11ex_;
    S_.i_syn_in_ *= V_.P11in_;

    // add evolution of presynaptic input current
    S_.i_syn_ex_ += ( 1. - V_.P11ex_ ) * S_.i_1_;

    // get read access to the correct input-buffer slot
    const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo( lag );
    auto& input = B_.input_buffer_.get_values_all_channels( input_buffer_slot );

    // the spikes arriving at T+1 have an immediate effect on the state of the neuron
    V_.weighted_spikes_ex_ = input[ Buffers_::SYN_EX ];
    V_.weighted_spikes_in_ = input[ Buffers_::SYN_IN ];

    S_.i_syn_ex_ += V_.weighted_spikes_ex_;
    S_.i_syn_in_ += V_.weighted_spikes_in_;

    double z = 0.0;  // spike state variable

    S_.surrogate_gradient_ =
      ( this->*compute_surrogate_gradient_ )( S_.r_ref_, S_.V_m_, P_.Theta_, P_.beta_, P_.gamma_ );

    if ( S_.V_m_ >= P_.Theta_ )
    {
      S_.r_ref_ = V_.RefractoryCounts_;
      S_.V_m_ = P_.V_reset_;

      set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

      SpikeEvent se;
      kernel().event_delivery_manager.send( *this, se, lag );

      z = 1.0;
    }

    append_new_eprop_history_entry( t );
    write_surrogate_gradient_to_history( t, S_.surrogate_gradient_ );
    write_firing_rate_reg_to_history( t, z, P_.f_target_, P_.kappa_reg_, P_.c_reg_ );

    S_.learning_signal_ = get_learning_signal_from_history( t );

    // set new input current
    S_.i_0_ = input[ Buffers_::I0 ];
    S_.i_1_ = input[ Buffers_::I1 ];

    // reset all values in the currently processed input-buffer slot
    B_.input_buffer_.reset_values_all_channels( input_buffer_slot );

    // log state data
    B_.logger_.record_data( t );
  }
}

/* ----------------------------------------------------------------
 * Event handling functions
 * ---------------------------------------------------------------- */

void
eprop_iaf_psc_exp::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ) );

  const double s = e.get_weight() * e.get_multiplicity();

  // separate buffer channels for excitatory and inhibitory inputs
  B_.input_buffer_.add_value( input_buffer_slot, s > 0 ? Buffers_::SYN_EX : Buffers_::SYN_IN, s );
}

void
eprop_iaf_psc_exp::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ) );

  if ( 0 == e.get_rport() )
  {
    B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I0, w * c );
  }
  if ( 1 == e.get_rport() )
  {
    B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I1, w * c );
  }
}

void
eprop_iaf_psc_exp::handle( LearningSignalConnectionEvent& e )
{
  for ( auto it_event = e.begin(); it_event != e.end(); )
  {
    const long time_step = e.get_stamp().get_steps();
    const double weight = e.get_weight();
    const double error_signal = e.get_coeffvalue( it_event );  // get_coeffvalue advances iterator
    const double learning_signal = weight * error_signal;

    write_learning_signal_to_history( time_step, learning_signal );
  }
}

void
eprop_iaf_psc_exp::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
eprop_iaf_psc_exp::compute_gradient( const long t_spike,
  const long t_spike_previous,
  double& z_previous_buffer,
  double& z_bar,
  double& e_bar,
  double& e_bar_reg,
  double& epsilon,
  double& weight,
  const CommonSynapseProperties& cp,
  WeightOptimizer* optimizer )
{
  double e = 0.0;                 // eligibility trace
  double z = 0.0;                 // spiking variable
  double z_current_buffer = 1.0;  // buffer containing the spike that triggered the current integration
  double psi = 0.0;               // surrogate gradient
  double L = 0.0;                 // learning signal
  double firing_rate_reg = 0.0;   // firing rate regularization
  double grad = 0.0;              // gradient

  const EpropSynapseCommonProperties& ecp = static_cast< const EpropSynapseCommonProperties& >( cp );
  const auto optimize_each_step = ( *ecp.optimizer_cp_ ).optimize_each_step_;

  auto eprop_hist_it = get_eprop_history( t_spike_previous - 1 );

  const long t_compute_until = std::min( t_spike_previous + V_.eprop_isi_trace_cutoff_steps_, t_spike );

  for ( long t = t_spike_previous; t < t_compute_until; ++t, ++eprop_hist_it )
  {
    z = z_previous_buffer;
    z_previous_buffer = z_current_buffer;
    z_current_buffer = 0.0;

    psi = eprop_hist_it->surrogate_gradient_;
    L = eprop_hist_it->learning_signal_;
    firing_rate_reg = eprop_hist_it->firing_rate_reg_;

    z_bar = V_.P22_ * z_bar + z;
    e = psi * z_bar;
    e_bar = P_.kappa_ * e_bar + e;
    e_bar_reg = P_.kappa_reg_ * e_bar_reg + ( 1.0 - P_.kappa_reg_ ) * e;

    if ( optimize_each_step )
    {
      grad = L * e_bar + firing_rate_reg * e_bar_reg;
      weight = optimizer->optimized_weight( *ecp.optimizer_cp_, t, grad, weight );
    }
    else
    {
      grad += L * e_bar + firing_rate_reg * e_bar_reg;
    }
  }

  if ( not optimize_each_step )
  {
    weight = optimizer->optimized_weight( *ecp.optimizer_cp_, t_compute_until, grad, weight );
  }

  const long cutoff_to_spike_interval = t_spike - t_compute_until;

  if ( cutoff_to_spike_interval > 0 )
  {
    z_bar *= std::pow( V_.P22_, cutoff_to_spike_interval );
    e_bar *= std::pow( P_.kappa_, cutoff_to_spike_interval );
    e_bar_reg *= std::pow( P_.kappa_reg_, cutoff_to_spike_interval );
  }
}

}  // namespace nest





// //// TAKEN FROM iaf_psc_exp
// /*
//  *  eprop_iaf_psc_exp.cpp
//  *
//  *  This file is part of NEST.
//  *
//  *  Copyright (C) 2004 The NEST Initiative
//  *
//  *  NEST is free software: you can redistribute it and/or modify
//  *  it under the terms of the GNU General Public License as published by
//  *  the Free Software Foundation, either version 2 of the License, or
//  *  (at your option) any later version.
//  *
//  *  NEST is distributed in the hope that it will be useful,
//  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  *  GNU General Public License for more details.
//  *
//  *  You should have received a copy of the GNU General Public License
//  *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
//  *
//  */

// #include "eprop_iaf_psc_exp.h"


// // Includes from libnestutil:
// #include "dict_util.h"
// #include "iaf_propagator.h"
// #include "numerics.h"

// // Includes from nestkernel:
// // #include "eprop_archiving_node_recurrent_impl.h"
// #include "exceptions.h"
// #include "iaf_propagator.h"
// #include "kernel_manager.h"
// #include "nest_impl.h"
// // #include "numerics.h"
// #include "ring_buffer_impl.h"
// #include "universal_data_logger_impl.h"

// /* ----------------------------------------------------------------
//  * Recordables map
//  * ---------------------------------------------------------------- */

// nest::RecordablesMap< nest::eprop_iaf_psc_exp > nest::eprop_iaf_psc_exp::recordablesMap_;

// namespace nest
// {
// void
// register_eprop_iaf_psc_exp( const std::string& name )
// {
//   register_node_model< eprop_iaf_psc_exp >( name );
// }

// // Override the create() method with one call to RecordablesMap::insert_()
// // for each quantity to be recorded.
// template <>
// void
// RecordablesMap< eprop_iaf_psc_exp >::create()
// {
//   // use standard names wherever you can for consistency!
//   insert_( names::V_m, &eprop_iaf_psc_exp::get_V_m_ );
//   insert_( names::I_syn_ex, &eprop_iaf_psc_exp::get_I_syn_ex_ );
//   insert_( names::I_syn_in, &eprop_iaf_psc_exp::get_I_syn_in_ );
// }
// }

// /* ----------------------------------------------------------------
//  * Default constructors defining default parameters and state
//  * ---------------------------------------------------------------- */

// nest::eprop_iaf_psc_exp::Parameters_::Parameters_()
//   : Tau_( 10.0 )              // in ms
//   , C_( 250.0 )               // in pF
//   , t_ref_( 2.0 )             // in ms
//   , E_L_( -70.0 )             // in mV
//   , I_e_( 0.0 )               // in pA
//   , Theta_( -55.0 - E_L_ )    // relative E_L_
//   , V_reset_( -70.0 - E_L_ )  // in mV
//   , tau_ex_( 2.0 )            // in ms
//   , tau_in_( 2.0 )            // in ms
//   , rho_( 0.2 )              // in 1/s
//   , delta_( 0.1 )             // in mV
// {
// }

// nest::eprop_iaf_psc_exp::State_::State_()
//   : i_0_( 0.0 )
//   , i_1_( 0.0 )
//   , i_syn_ex_( 0.0 )
//   , i_syn_in_( 0.0 )
//   , V_m_( 0.0 )
//   , r_ref_( 0 )
// {
// }

// /* ----------------------------------------------------------------
//  * Parameter and state extractions and manipulation functions
//  * ---------------------------------------------------------------- */

// void
// nest::eprop_iaf_psc_exp::Parameters_::get( Dictionary& d ) const
// {
//   d[ names::E_L ] = E_L_;  // resting potential
//   d[ names::I_e ] = I_e_;
//   d[ names::V_th ] = Theta_ + E_L_;  // threshold value
//   d[ names::V_reset ] = V_reset_ + E_L_;
//   d[ names::C_m ] = C_;
//   d[ names::tau_m ] = Tau_;
//   d[ names::tau_syn_ex ] = tau_ex_;
//   d[ names::tau_syn_in ] = tau_in_;
//   d[ names::t_ref ] = t_ref_;
//   d[ names::rho ] = rho_;
//   d[ names::delta ] = delta_;
// }

// double
// nest::eprop_iaf_psc_exp::Parameters_::set( const Dictionary& d, Node* node )
// {
//   // if E_L_ is changed, we need to adjust all variables defined relative to
//   // E_L_
//   const double ELold = E_L_;
//   update_value_param( d, names::E_L, E_L_, node );
//   const double delta_EL = E_L_ - ELold;

//   if ( update_value_param( d, names::V_reset, V_reset_, node ) )
//   {
//     V_reset_ -= E_L_;
//   }
//   else
//   {
//     V_reset_ -= delta_EL;
//   }

//   if ( update_value_param( d, names::V_th, Theta_, node ) )
//   {
//     Theta_ -= E_L_;
//   }
//   else
//   {
//     Theta_ -= delta_EL;
//   }

//   update_value_param( d, names::I_e, I_e_, node );
//   update_value_param( d, names::C_m, C_, node );
//   update_value_param( d, names::tau_m, Tau_, node );
//   update_value_param( d, names::tau_syn_ex, tau_ex_, node );
//   update_value_param( d, names::tau_syn_in, tau_in_, node );
//   update_value_param( d, names::t_ref, t_ref_, node );
//   if ( V_reset_ >= Theta_ )
//   {
//     throw BadProperty( "Reset potential must be smaller than threshold." );
//   }
//   if ( C_ <= 0 )
//   {
//     throw BadProperty( "Capacitance must be strictly positive." );
//   }
//   if ( Tau_ <= 0 or tau_ex_ <= 0 or tau_in_ <= 0 )
//   {
//     throw BadProperty( "Membrane and synapse time constants must be strictly positive." );
//   }
//   if ( t_ref_ < 0 )
//   {
//     throw BadProperty( "Refractory time must not be negative." );
//   }

//   d.update_value( "rho", rho_ );
//   if ( rho_ < 0 )
//   {
//     throw BadProperty( "Stochastic firing intensity must not be negative." );
//   }

//   d.update_value( "delta", delta_ );
//   if ( delta_ < 0 )
//   {
//     throw BadProperty( "Width of threshold region must not be negative." );
//   }

//   return delta_EL;
// }

// void
// nest::eprop_iaf_psc_exp::State_::get( Dictionary& d, const Parameters_& p ) const
// {
//   d[ names::V_m ] = V_m_ + p.E_L_;  // Membrane potential
// }

// void
// nest::eprop_iaf_psc_exp::State_::set( const Dictionary& d, const Parameters_& p, double delta_EL, Node* node )
// {
//   if ( update_value_param( d, names::V_m, V_m_, node ) )
//   {
//     V_m_ -= p.E_L_;
//   }
//   else
//   {
//     V_m_ -= delta_EL;
//   }
// }

// nest::eprop_iaf_psc_exp::Buffers_::Buffers_( eprop_iaf_psc_exp& n )
//   : logger_( n )
// {
// }

// nest::eprop_iaf_psc_exp::Buffers_::Buffers_( const Buffers_&, eprop_iaf_psc_exp& n )
//   : logger_( n )
// {
// }

// /* ----------------------------------------------------------------
//  * Default and copy constructor for node
//  * ---------------------------------------------------------------- */

// nest::eprop_iaf_psc_exp::eprop_iaf_psc_exp()
//   : ArchivingNode()
//   , P_()
//   , S_()
//   , B_( *this )
// {
//   recordablesMap_.create();
// }

// nest::eprop_iaf_psc_exp::eprop_iaf_psc_exp( const eprop_iaf_psc_exp& n )
//   : ArchivingNode( n )
//   , P_( n.P_ )
//   , S_( n.S_ )
//   , B_( n.B_, *this )
// {
// }

// /* ----------------------------------------------------------------
//  * Node initialization functions
//  * ---------------------------------------------------------------- */

// void
// nest::eprop_iaf_psc_exp::init_buffers_()
// {
//   B_.input_buffer_.clear();  // includes resize
//   B_.logger_.reset();
//   ArchivingNode::clear_history();
// }

// void
// nest::eprop_iaf_psc_exp::pre_run_hook()
// {
//   // ensures initialization in case mm connected after Simulate
//   B_.logger_.init();

//   const double h = Time::get_resolution().get_ms();

//   // these P are independent
//   V_.P11ex_ = std::exp( -h / P_.tau_ex_ );
//   V_.P11in_ = std::exp( -h / P_.tau_in_ );

//   V_.P22_ = std::exp( -h / P_.Tau_ );

//   // these are determined according to a numeric stability criterion
//   V_.P21ex_ = IAFPropagatorExp( P_.tau_ex_, P_.Tau_, P_.C_ ).evaluate( h );
//   V_.P21in_ = IAFPropagatorExp( P_.tau_in_, P_.Tau_, P_.C_ ).evaluate( h );

//   V_.P20_ = P_.Tau_ / P_.C_ * ( 1.0 - V_.P22_ );

//   // t_ref_ specifies the length of the absolute refractory period as
//   // a double in ms. The grid based eprop_iaf_psc_exp can only handle refractory
//   // periods that are integer multiples of the computation step size (h).
//   // To ensure consistency with the overall simulation scheme such conversion
//   // should be carried out via objects of class nest::Time. The conversion
//   // requires 2 steps:
//   //     1. A time object r is constructed, defining representation of
//   //        t_ref_ in tics. This representation is then converted to computation
//   //        time steps again by a strategy defined by class nest::Time.
//   //     2. The refractory time in units of steps is read out get_steps(), a
//   //        member function of class nest::Time.
//   //
//   // Choosing a t_ref_ that is not an integer multiple of the computation time
//   // step h will lead to accurate (up to the resolution h) and self-consistent
//   // results. However, a neuron model capable of operating with real valued
//   // spike time may exhibit a different effective refractory time.

//   V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
//   // since t_ref_ >= 0, this can only fail in error
//   assert( V_.RefractoryCounts_ >= 0 );

//   V_.rng_ = get_vp_specific_rng( get_thread() );
// }

// void
// nest::eprop_iaf_psc_exp::update( const Time& origin, const long from, const long to )
// {
//   const double h = Time::get_resolution().get_ms();

//   // evolve from timestep 'from' to timestep 'to' with steps of h each
//   for ( long lag = from; lag < to; ++lag )
//   {
//     if ( S_.r_ref_ == 0 )  // neuron not refractory, so evolve V
//     {
//       S_.V_m_ =
//         S_.V_m_ * V_.P22_ + S_.i_syn_ex_ * V_.P21ex_ + S_.i_syn_in_ * V_.P21in_ + ( P_.I_e_ + S_.i_0_ ) * V_.P20_;
//     }
//     else
//     {
//       // neuron is absolute refractory
//       --S_.r_ref_;
//     }

//     // exponential decaying PSCs
//     S_.i_syn_ex_ *= V_.P11ex_;
//     S_.i_syn_in_ *= V_.P11in_;

//     // add evolution of presynaptic input current
//     S_.i_syn_ex_ += ( 1. - V_.P11ex_ ) * S_.i_1_;

//     // get read access to the correct input-buffer slot
//     const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo( lag );
//     auto& input = B_.input_buffer_.get_values_all_channels( input_buffer_slot );

//     // the spikes arriving at T+1 have an immediate effect on the state of the
//     // neuron

//     V_.weighted_spikes_ex_ = input[ Buffers_::SYN_EX ];
//     V_.weighted_spikes_in_ = input[ Buffers_::SYN_IN ];

//     S_.i_syn_ex_ += V_.weighted_spikes_ex_;
//     S_.i_syn_in_ += V_.weighted_spikes_in_;

//     if ( ( P_.delta_ < 1e-10 and S_.V_m_ >= P_.Theta_ )                    // deterministic threshold crossing
//       or ( P_.delta_ > 1e-10 and V_.rng_->drand() < phi_() * h * 1e-3 ) )  // stochastic threshold crossing
//     {
//       S_.r_ref_ = V_.RefractoryCounts_;
//       S_.V_m_ = P_.V_reset_;

//       set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

//       SpikeEvent se;
//       kernel().event_delivery_manager.send( *this, se, lag );
//     }

//     // set new input current
//     S_.i_0_ = input[ Buffers_::I0 ];
//     S_.i_1_ = input[ Buffers_::I1 ];

//     // reset all values in the currently processed input-buffer slot
//     B_.input_buffer_.reset_values_all_channels( input_buffer_slot );

//     // log state data
//     B_.logger_.record_data( origin.get_steps() + lag );
//   }
// }

// void
// nest::eprop_iaf_psc_exp::handle( SpikeEvent& e )
// {
//   assert( e.get_delay_steps() > 0 );

//   const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo(
//     e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ) );

//   const double s = e.get_weight() * e.get_multiplicity();

//   // separate buffer channels for excitatory and inhibitory inputs
//   B_.input_buffer_.add_value( input_buffer_slot, s > 0 ? Buffers_::SYN_EX : Buffers_::SYN_IN, s );
// }

// void
// nest::eprop_iaf_psc_exp::handle( CurrentEvent& e )
// {
//   assert( e.get_delay_steps() > 0 );

//   const double c = e.get_current();
//   const double w = e.get_weight();

//   const size_t input_buffer_slot = kernel().event_delivery_manager.get_modulo(
//     e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ) );

//   if ( 0 == e.get_rport() )
//   {
//     B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I0, w * c );
//   }
//   if ( 1 == e.get_rport() )
//   {
//     B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I1, w * c );
//   }
// }

// void
// nest::eprop_iaf_psc_exp::handle( DataLoggingRequest& e )
// {
//   B_.logger_.handle( e );
// }
