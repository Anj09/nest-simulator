#ifndef EPROP_IAF_PSC_EXP_H
#define EPROP_IAF_PSC_EXP_H

// nestkernel
#include "connection.h"
#include "eprop_archiving_node_impl.h"
#include "eprop_archiving_node_recurrent.h"
#include "eprop_synapse.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{

/*
##################      Neuron parameters      ##################

Parameter              Unit    Math equivalent         Default  
====================== ======= ======================= ======== 
C_m                    pF      :math:`C_m`             250.0    
E_L                    mV      :math:`E_L`             -70.0    
I_e                    pA      :math:`I_e`             0.0      
t_ref                  ms      :math:`t_{ref}`         2.0      
tau_m                  ms      :math:`\tau_m`          10.0     
tau_syn_ex             ms      :math:`\tau_{syn,ex}`   2.0      
tau_syn_in             ms      :math:`\tau_{syn,in}`   2.0      
V_th                   mV      :math:`v_{th}`          -55.0    
V_reset                mV      :math:`v_{reset}`       -70.0    



#################      E-prop parameters       #################

Parameter              Unit    Math equivalent            Default          
====================== ======= ========================== ================ 
beta                           :math:`\beta`              1.0              
c_reg                          :math:`c_{reg}`            0.0              
eprop_isi_trace_cutoff ms      :math:`\Delta t_c`         1000.0           
f_target               Hz      :math:`f^{target}`         10.0             
gamma                          :math:`\gamma`             0.3              
kappa                          :math:`\kappa`             0.97             
kappa_reg                      :math:`\kappa_{reg}`       0.97             
surrogate_gradient                                        piecewise_linear 


#########        Recordables       ###########

===================================================================
V_m               Membrane potential
I_syn_ex          Excitatory synaptic current
I_syn_in          Inhibitory synaptic current
learning_signal   Learning signal from readout neurons
surrogate_gradient Surrogate gradient (pseudo-derivative)
 */

void register_eprop_iaf_psc_exp( const std::string& name );


class eprop_iaf_psc_exp : public EpropArchivingNodeRecurrent< false >
{

public:
  //Default constructor.
  eprop_iaf_psc_exp();

  //Copy constructor.
  eprop_iaf_psc_exp( const eprop_iaf_psc_exp& );

  using Node::handle;
  using Node::handles_test_event;

  size_t send_test_event( Node&, size_t, synindex, bool ) override;

  void handle( SpikeEvent& ) override;
  void handle( CurrentEvent& ) override;
  void handle( LearningSignalConnectionEvent& ) override;
  void handle( DataLoggingRequest& ) override;

  size_t handles_test_event( SpikeEvent&, size_t ) override;
  size_t handles_test_event( CurrentEvent&, size_t ) override;
  size_t handles_test_event( LearningSignalConnectionEvent&, size_t ) override;
  size_t handles_test_event( DataLoggingRequest&, size_t ) override;

  void get_status( Dictionary& ) const override;
  void set_status( const Dictionary& ) override;

private:
  void init_buffers_() override;
  void pre_run_hook() override;

  void update( const Time&, const long, const long ) override;

  //Computing gradient for e-prop weight update 
  void compute_gradient( const long,
    const long,
    double&,
    double&,
    double&,
    double&,
    double&,
    double&,
    const CommonSynapseProperties&,
    WeightOptimizer*,
    const bool,
    const bool,
    double&,
    long&,
    long& ) override;
  // void compute_gradient( const long,
  //   const long,
  //   double&,
  //   double&,
  //   double&,
  //   double&,
  //   double&,
  //   double&,
  //   const CommonSynapseProperties&,
  //   WeightOptimizer* ) override;

  long get_shift() const override;
  bool is_eprop_recurrent_node() const override;
  // long get_eprop_isi_trace_cutoff() const override;

 
  friend class RecordablesMap< eprop_iaf_psc_exp >;


  friend class UniversalDataLogger< eprop_iaf_psc_exp >;

  // Parameters
 

  struct Parameters_
  {
    
    double Tau_;

    
    double C_;

   
    double t_ref_;

   
    double E_L_;

  
    double I_e_;

    
    double Theta_;

  
    double V_reset_;

    
    double tau_ex_;

    
    double tau_in_;

 
    double c_reg_;

    
    double f_target_;


    double beta_;


    double gamma_;

    //surogate gradient function 
    std::string surrogate_gradient_function_;

    
    double kappa_;

    
    double kappa_reg_;

    
    double eprop_isi_trace_cutoff_;

    Parameters_();

    void get( Dictionary& ) const;
    double set( const Dictionary&, Node* );
  };


  // State variables
  

  struct State_
  {
    
    double i_0_;

    
    double i_1_;

    
    double i_syn_ex_;

    
    double i_syn_in_;

    
    double V_m_;

    
    int r_ref_;

    
    double learning_signal_;

    
    double surrogate_gradient_;

    State_();

    void get( Dictionary&, const Parameters_& ) const;
    void set( const Dictionary&, const Parameters_&, double, Node* );
  };

  
  // Buffers

  struct Buffers_
  {
    Buffers_( eprop_iaf_psc_exp& );
    Buffers_( const Buffers_&, eprop_iaf_psc_exp& );

    
    enum
    {
      SYN_IN = 0,
      SYN_EX,
      I0,
      I1,
      NUM_INPUT_CHANNELS
    };

    
    MultiChannelInputBuffer< NUM_INPUT_CHANNELS > input_buffer_;

   
    UniversalDataLogger< eprop_iaf_psc_exp > logger_;
  };

  
  // Internal variables

  struct Variables_
  {

    double P22_;

   
    double P20_;

    
    double P11ex_;

    
    double P11in_;

    
    double P21ex_;

  
    double P21in_;

    
    double weighted_spikes_ex_;

   
    double weighted_spikes_in_;

   
    int RefractoryCounts_;


    long eprop_isi_trace_cutoff_steps_;
  };


  // Getter functions for recordables

  double
  get_V_m_() const
  {
    return S_.V_m_ + P_.E_L_;
  }

  double
  get_I_syn_ex_() const
  {
    return S_.i_syn_ex_;
  }

  double
  get_I_syn_in_() const
  {
    return S_.i_syn_in_;
  }

  double
  get_surrogate_gradient_() const
  {
    return S_.surrogate_gradient_;
  }

  double
  get_learning_signal_() const
  {
    return S_.learning_signal_;
  }

  // Member data

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  static RecordablesMap< eprop_iaf_psc_exp > recordablesMap_;
};


// Inline function

inline long
eprop_iaf_psc_exp::get_shift() const
{
  return offset_gen_ + delay_in_rec_;
}

inline bool
eprop_iaf_psc_exp::is_eprop_recurrent_node() const
{
  return true;
}

// inline long
// eprop_iaf_psc_exp::get_eprop_isi_trace_cutoff() const
// {
//   return V_.eprop_isi_trace_cutoff_steps_;
// }

inline size_t
eprop_iaf_psc_exp::send_test_event( Node& target, size_t receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline size_t
eprop_iaf_psc_exp::handles_test_event( SpikeEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline size_t
eprop_iaf_psc_exp::handles_test_event( CurrentEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline size_t
eprop_iaf_psc_exp::handles_test_event( LearningSignalConnectionEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline size_t
eprop_iaf_psc_exp::handles_test_event( DataLoggingRequest& dlr, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
eprop_iaf_psc_exp::get_status( Dictionary& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  d[ names::recordables ] = recordablesMap_.get_list();
}

inline void
eprop_iaf_psc_exp::set_status( const Dictionary& d )
{
  Parameters_ ptmp = P_;
  State_ stmp = S_;

  const double delta_EL = ptmp.set( d, this );
  stmp.set( d, ptmp, delta_EL, this );

  P_ = ptmp;
  S_ = stmp;
}

}  // namespace nest

#endif  // EPROP_IAF_PSC_EXP_H
