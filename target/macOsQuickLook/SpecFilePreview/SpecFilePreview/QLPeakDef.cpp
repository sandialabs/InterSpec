#include <memory>
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/poisson.hpp>


#include "SpecUtils/3rdparty/rapidxml/rapidxml.hpp"
#include "SpecUtils/3rdparty/rapidxml/rapidxml_utils.hpp"
#include "SpecUtils/3rdparty/rapidxml/rapidxml_print.hpp"

   
#include "QLPeakDef.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"



using namespace std;

#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH


const int QLPeakDef::sm_xmlSerializationVersion = 0;
const int QLPeakContinuum::sm_xmlSerializationVersion = 0;


namespace
{
  //clones 'source' into the document that 'result' is a part of.
  //  'result' is cleared and set lexically equal to 'source'.
  void clone_node_deep( const ::rapidxml::xml_node<char> *source,
                        ::rapidxml::xml_node<char> *result )
  {
    using namespace ::rapidxml;
    
    xml_document<char> *doc = result->document();
    if( !doc )
      throw runtime_error( "clone_node_deep: insert result into document before calling" );
    
    result->remove_all_attributes();
    result->remove_all_nodes();
    result->type(source->type());
    
    
    // Clone name and value
    char *str = doc->allocate_string( source->name(), source->name_size() );
    result->name( str, source->name_size());
    
    if( source->value() )
    {
      str = doc->allocate_string( source->value(), source->value_size() );
      result->value( str, source->value_size() );
    }
    
    // Clone child nodes and attributes
    for( xml_node<char> *child = source->first_node(); child; child = child->next_sibling() )
    {
      xml_node<char> *clone = doc->allocate_node( child->type() );
      result->append_node( clone );
      clone_node_deep( child, clone );
    }
    
    for( xml_attribute<char> *attr = source->first_attribute(); attr; attr = attr->next_attribute())
    {
      const char *name = doc->allocate_string( attr->name(), attr->name_size() );
      const char *value = doc->allocate_string( attr->value(), attr->value_size() );
      xml_attribute<char> *clone = doc->allocate_attribute(name, value, attr->name_size(), attr->value_size());
      result->append_attribute( clone );
    }
  }//void clone_node_deep(...)
  
  double landau_cdf(double x, double xi, double x0) {
    // implementation of landau distribution (from DISLAN)
    //The algorithm was taken from the Cernlib function dislan(G110)
    //Reference: K.S.Kolbig and B.Schorr, "A program package for the Landau
    //distribution", Computer Phys.Comm., 31(1984), 97-111
    //
    //Lifted from the root/math/mathcore/src/ProbFuncMathCore.cxx file
    //  by wcjohns 20120216
    
    static double p1[5] = {0.2514091491e+0,-0.6250580444e-1, 0.1458381230e-1, -0.2108817737e-2, 0.7411247290e-3};
    static double q1[5] = {1.0            ,-0.5571175625e-2, 0.6225310236e-1, -0.3137378427e-2, 0.1931496439e-2};
    
    static double p2[4] = {0.2868328584e+0, 0.3564363231e+0, 0.1523518695e+0, 0.2251304883e-1};
    static double q2[4] = {1.0            , 0.6191136137e+0, 0.1720721448e+0, 0.2278594771e-1};
    
    static double p3[4] = {0.2868329066e+0, 0.3003828436e+0, 0.9950951941e-1, 0.8733827185e-2};
    static double q3[4] = {1.0            , 0.4237190502e+0, 0.1095631512e+0, 0.8693851567e-2};
    
    static double p4[4] = {0.1000351630e+1, 0.4503592498e+1, 0.1085883880e+2, 0.7536052269e+1};
    static double q4[4] = {1.0            , 0.5539969678e+1, 0.1933581111e+2, 0.2721321508e+2};
    
    static double p5[4] = {0.1000006517e+1, 0.4909414111e+2, 0.8505544753e+2, 0.1532153455e+3};
    static double q5[4] = {1.0            , 0.5009928881e+2, 0.1399819104e+3, 0.4200002909e+3};
    
    static double p6[4] = {0.1000000983e+1, 0.1329868456e+3, 0.9162149244e+3, -0.9605054274e+3};
    static double q6[4] = {1.0            , 0.1339887843e+3, 0.1055990413e+4, 0.5532224619e+3};
    
    static double a1[4] = {0, -0.4583333333e+0, 0.6675347222e+0,-0.1641741416e+1};
    
    static double a2[4] = {0,  1.0            ,-0.4227843351e+0,-0.2043403138e+1};
    
    double v = (x - x0)/xi;
    double u;
    double lan;
    
    if (v < -5.5) {
      u = std::exp(v+1);
      lan = 0.3989422803*std::exp(-1./u)*std::sqrt(u)*(1+(a1[1]+(a1[2]+a1[3]*u)*u)*u);
    }
    else if (v < -1 ) {
      u = std::exp(-v-1);
      lan = (std::exp(-u)/std::sqrt(u))*(p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v)*v)*v)*v)/
      (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v)*v)*v)*v);
    }
    else if (v < 1)
      lan = (p2[0]+(p2[1]+(p2[2]+p2[3]*v)*v)*v)/(q2[0]+(q2[1]+(q2[2]+q2[3]*v)*v)*v);
    else if (v < 4)
      lan = (p3[0]+(p3[1]+(p3[2]+p3[3]*v)*v)*v)/(q3[0]+(q3[1]+(q3[2]+q3[3]*v)*v)*v);
    else if (v < 12) {
      u = 1./v;
      lan = (p4[0]+(p4[1]+(p4[2]+p4[3]*u)*u)*u)/(q4[0]+(q4[1]+(q4[2]+q4[3]*u)*u)*u);
    }
    else if (v < 50) {
      u = 1./v;
      lan = (p5[0]+(p5[1]+(p5[2]+p5[3]*u)*u)*u)/(q5[0]+(q5[1]+(q5[2]+q5[3]*u)*u)*u);
    }
    else if (v < 300) {
      u = 1./v;
      lan = (p6[0]+(p6[1]+(p6[2]+p6[3]*u)*u)*u)/(q6[0]+(q6[1]+(q6[2]+q6[3]*u)*u)*u);
    }
    else {
      u = 1./(v-v*std::log(v)/(v+1));
      lan = 1-(a2[1]+(a2[2]+a2[3]*u)*u)*u;
    }
    
    return lan;
  }//double landau_cdf(double x, double xi, double x0)

}//namespace


double skewedGaussianIntegral( double x0, double x1,
                               double mu, double s,
                               double L )
{
  using boost::math::erf;
  using boost::math::erfc;
  static const double sqrt2 = boost::math::constants::root_two<double>();
  
  return -0.5*erf((x0 - mu)/(sqrt2*s)) + erf((x1 - mu)/(sqrt2*s))
   + exp((L*(L*s*s - 2*x0 + 2*mu))/2)*erfc((L*s*s - x0 + mu)/(sqrt2*s))
  - exp((L*(L*s*s - 2*x1 + 2*mu))/2) *erfc((L*s*s - x1 + mu)/(sqrt2*s));
}//double skewedGaussianIntegral(...)

double skewedGaussianIndefinitIntegral( double x,
                              double A, double c,
                              double w, double t )
{
  using boost::math::erf;
  using boost::math::erfc;
  static const double sqrt2 = boost::math::constants::root_two<double>();
  
//  integral (A e^(1/2 (w/t)^2-(x-c)/t) (1/2+1/2 erf(((x-c)/w-w/t)/sqrt(2))))/t dx =
  return -0.5*A*exp(-x/t)*(exp((2*c*t+w*w)/(2*t*t))*(erf(((x-c)/w-w/t)/sqrt2)+1)-exp(x/t)*erf((x-c)/(sqrt2*w)));
}

double skewedGaussianIntegral( double x0, double x1,
                              double A, double c,
                              double w, double t )
{
  return skewedGaussianIndefinitIntegral( x1, A, c, w, t ) - skewedGaussianIndefinitIntegral( x0, A, c, w, t );
}



QLPeakDef::QLPeakDef()
{
  reset();
}


void QLPeakDef::reset()
{
  m_userLabel                = "";
  m_type                     = GaussianDefined;
  m_skewType                 = QLPeakDef::NoSkew;

  m_parentNuclide            = "";
  m_gammaEnergy              = 0.0;
  //m_transition               = NULL;
  //m_radparticleIndex         = -1;
  m_sourceGammaType          = NormalGamma;
  m_useForCalibration        = true;
  m_useForShieldingSourceFit = false;
  
  m_xrayElement              = "";
  m_xrayEnergy               = 0.0;
  m_reaction                 = "";
  m_reactionEnergy           = 0.0;

  std::shared_ptr<QLPeakContinuum> newcont = std::make_shared<QLPeakContinuum>();
  m_continuum = newcont;
  
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    m_coefficients[t] = 0.0;
    m_uncertainties[t] = -1.0;
    
    switch( t )
    {
      case QLPeakDef::Mean:
      case QLPeakDef::Sigma:
      case QLPeakDef::GaussAmplitude:
        m_fitFor[t] = true;
      break;
        
      case QLPeakDef::LandauAmplitude:
      case QLPeakDef::LandauMode:
      case QLPeakDef::LandauSigma:
      case QLPeakDef::Chi2DOF:
      case QLPeakDef::NumCoefficientTypes:
        m_fitFor[t] = false;
      break;
    }//switch( type )
  }//for( loop over coefficients )
}//void QLPeakDef::reset()


QLPeakDef::QLPeakDef( double m, double s, double a )
{
  reset();
  m_coefficients[QLPeakDef::Mean] = m;
  m_coefficients[QLPeakDef::Sigma] = s;
  m_coefficients[QLPeakDef::GaussAmplitude] = a;
}



QLPeakDef::QLPeakDef( double xlow, double xhigh, double mean,
                    std::shared_ptr<const Measurement> data, std::shared_ptr<const Measurement> background )
{
  reset();
  m_type = QLPeakDef::DataDefined;
  m_coefficients[QLPeakDef::Mean] = mean;
  m_continuum->setRange( xlow, xhigh );
  
  if( !data )
    return;

  m_continuum->setType( QLPeakContinuum::External );
  m_continuum->setExternalContinuum( background );
  m_continuum->setRange( xlow, xhigh );
  
  m_coefficients[QLPeakDef::GaussAmplitude] = gamma_integral( data, xlow, xhigh );
  if( background )
    m_coefficients[QLPeakDef::GaussAmplitude] -= gamma_integral( background, xlow, xhigh );
}//QLPeakDef( constructor )



const char *QLPeakDef::to_string( const CoefficientType type )
{
  switch( type )
  {
    case QLPeakDef::Mean:                return "Centroid";
    case QLPeakDef::Sigma:               return "Width";
    case QLPeakDef::GaussAmplitude:      return "Amplitude";
    case QLPeakDef::LandauAmplitude:     return "LandauAmplitude";
    case QLPeakDef::LandauMode:          return "LandauMode";
    case QLPeakDef::LandauSigma:         return "LandauSigma";
    case QLPeakDef::Chi2DOF:             return "Chi2";
    case QLPeakDef::NumCoefficientTypes: return "";
  }//switch( type )

  return "";
}//const char *QLPeakDef::to_string( const CoefficientType type )


void QLPeakContinuum::toXml( rapidxml::xml_node<char> *parent, const int contId ) const
{
  using namespace rapidxml;
  
  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  
  if( !doc )
    throw runtime_error( "QLPeakContinuum::toXml(...): invalid input" );
  
  char buffer[128];
  xml_node<char> *node = 0;
  xml_node<char> *cont_node = doc->allocate_node( node_element, "PeakContinuum" );

  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  cont_node->append_attribute( att );

  snprintf( buffer, sizeof(buffer), "%i", contId );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "id", val );
  cont_node->append_attribute( att );
    
  parent->append_node( cont_node );
  
  const char *type = 0;
  switch( m_type )
  {
    case NoOffset:   type = "NoOffset";   break;
    case Constant:   type = "Constant";   break;
    case Linear:     type = "Linear";     break;
    case Quardratic: type = "Quardratic"; break;
    case Cubic:      type = "Cubic";      break;
    case External:   type = "External";   break;
  }//switch( m_type )
  
  node = doc->allocate_node( node_element, "Type", type );
  cont_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8e", m_lowerEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "LowerEnergy", val );
  cont_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8e", m_upperEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "UpperEnergy", val );
  cont_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8e", m_refernceEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "ReferenceEnergy", val );
  cont_node->append_node( node );
  
  if( m_type != NoOffset && m_type != External )
  {
    stringstream valsstrm, uncertstrm, fitstrm;
    for( size_t i = 0; i < m_values.size(); ++i )
    {
      const char *spacer = (i ? " " : "");
      snprintf( buffer, sizeof(buffer), "%1.8e", m_values[i] );  
      valsstrm << spacer << buffer;
    
      snprintf( buffer, sizeof(buffer), "%1.8e", m_uncertainties[i] );  
      uncertstrm << spacer << buffer;
    
      fitstrm << spacer << (m_fitForValue[i] ? '1': '0');
    }//for( size_t i = 0; i < m_values.size(); ++i )
    
    xml_node<char> *coeffs_node = doc->allocate_node( node_element, "Coefficients" );
    cont_node->append_node( coeffs_node );
        
    val = doc->allocate_string( valsstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Values", val );
    coeffs_node->append_node( node );
    
    val = doc->allocate_string( uncertstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Uncertainties", val );
    coeffs_node->append_node( node );
    
    val = doc->allocate_string( fitstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Fittable", val );
    coeffs_node->append_node( node );
  }//if( m_type != NoOffset && m_type != External )
  
  if( !!m_externalContinuum )
  {
    stringstream contXml;
    m_externalContinuum->write_2006_N42_xml( contXml );
    //We actually need to parse the XML here, and then insert it into the hierarchy
    
    const string datastr = contXml.str();
    std::unique_ptr<char[]> data( new char [datastr.size()+1] );
    //boost::scoped_array<char> data( new char [datastr.size()+1] );
    memcpy( data.get(), datastr.c_str(), datastr.size()+1 );
    
    xml_document<char> contdoc;
    const int flags = rapidxml::parse_normalize_whitespace
                     | rapidxml::parse_trim_whitespace;
    contdoc.parse<flags>( data.get() );
    
    node = doc->allocate_node( node_element, "ExternalContinuum", val );
    cont_node->append_node( node );
    
    xml_node<char> *spec_node = contdoc.first_node( "Measurement", 11 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Measurement node" );
    spec_node = spec_node->first_node( "Spectrum", 8 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Spectrum node" );
    
    xml_node<char> *new_spec_node = doc->allocate_node( node_element );
    node->append_node( new_spec_node );
    
    clone_node_deep( spec_node, new_spec_node );
  }//if( !!m_externalContinuum )
}//void QLPeakContinuum::toXml(...)





void QLPeakContinuum::fromXml( const rapidxml::xml_node<char> *cont_node, int &contId )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !cont_node )
    throw runtime_error( "PeakContinuum::fromXml(...): invalid input" );
  
  if( !compare( cont_node->name(), cont_node->name_size(), "PeakContinuum", 13, false ) )
    throw std::logic_error( "QLPeakContinuum::fromXml(...): invalid input node name" );
  
  xml_attribute<char> *att = cont_node->first_attribute( "version", 7 );
  
  int version;
  if( !att || !att->value() || (sscanf(att->value(), "%i", &version)!=1) )
    throw runtime_error( "QLPeakContinuum invalid version" );
  
  if( version != sm_xmlSerializationVersion )
    throw runtime_error( "Invalid QLPeakContinuum version" );
  
  att = cont_node->first_attribute( "id", 2 );
  if( !att || !att->value() || (sscanf(att->value(), "%i", &contId)!=1) )
    throw runtime_error( "QLPeakContinuum invalid ID" );

  xml_node<char> *node = cont_node->first_node( "Type", 4 );

  if( !node || !node->value() )
    throw runtime_error( "QLPeakContinuum not Type node" );
  
  if( compare(node->value(),node->value_size(),"NoOffset",8,false) )
    m_type = NoOffset;
  else if( compare(node->value(),node->value_size(),"Constant",8,false) )
    m_type = Constant;
  else if( compare(node->value(),node->value_size(),"Linear",6,false) )
    m_type = Linear;
  else if( compare(node->value(),node->value_size(),"Quardratic",10,false) )
    m_type = Quardratic;
  else if( compare(node->value(),node->value_size(),"Cubic",5,false) )
    m_type = Cubic;
  else if( compare(node->value(),node->value_size(),"External",8,false) )
    m_type = External;
  else
    throw runtime_error( "Invalid continuum type" );
  
  float dummyval;
  node = cont_node->first_node( "LowerEnergy", 11 );
  if( !node || !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
    throw runtime_error( "Continuum didnt have LowerEnergy" );
  m_lowerEnergy = dummyval;
    
  node = cont_node->first_node( "UpperEnergy", 11 );
  if( !node || !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
    throw runtime_error( "Continuum didnt have UpperEnergy" );
  m_upperEnergy = dummyval;
  
  node = cont_node->first_node( "ReferenceEnergy", 15 );
  if( !node || !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
    throw runtime_error( "Continuum didnt have ReferenceEnergy" );
  m_refernceEnergy = dummyval;
  
  if( m_type != NoOffset && m_type != External )
  {
    xml_node<char> *coeffs_node = cont_node->first_node("Coefficients",12);
    if( !coeffs_node )
      throw runtime_error( "Contoinuum didt have Coefficients node" );
    
    std::vector<float> contents;
    node = coeffs_node->first_node( "Values", 6 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Values" );
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_values.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_values[i] = contents[i]; 
    
    
    node = coeffs_node->first_node( "Uncertainties", 13 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Uncertainties" );  
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_uncertainties.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_uncertainties[i] = contents[i]; 
    
    
    node = coeffs_node->first_node( "Fittable", 8 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Fittable" );  
    
    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_fitForValue.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_fitForValue[i] = (contents[i] > 0.5f); 
    
    if( m_values.size() != m_uncertainties.size() 
        || m_fitForValue.size() != m_values.size() )
      throw runtime_error( "Continuum coefficients not consistent" );
  }else
  {
    m_values.clear();
    m_uncertainties.clear();
    m_fitForValue.clear();
  }//if( m_type != NoOffset && m_type != External ) / else
  
  
  node = cont_node->first_node( "ExternalContinuum", 17 );
  if( node )
  {
    node = node->first_node( "Spectrum", 8 );
    if( !node )
      throw runtime_error( "Spectrum node expected under ExternalContinuum" );
    std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
    m_externalContinuum = meas;
    meas->set_2006_N42_spectrum_node_info( node );
    
    if( !meas->channel_energies() || meas->channel_energies()->empty() )
      meas->popuplate_channel_energies_from_coeffs();    
  }//if( node )
}//void QLPeakContinuum::fromXml(...)



rapidxml::xml_node<char> *QLPeakDef::toXml( rapidxml::xml_node<char> *parent,
                     rapidxml::xml_node<char> *continuum_parent,
           std::map<std::shared_ptr<QLPeakContinuum>,int> &continuums ) const
{
  using namespace rapidxml;
  
  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  
  if( !doc )
    throw runtime_error( "QLPeakDef::toXml(...): invalid input" );
  
  if( !m_continuum )
    throw logic_error( "QLPeakDef::toXml(...): continuum should be valid" );
  
  if( !continuums.count(m_continuum) )
  {
    const int index = static_cast<int>( continuums.size() + 1 );
    m_continuum->toXml( continuum_parent, index );
    continuums[m_continuum] = index;
  }//if( !continuums.count(m_continuum) )
  
  char buffer[128];
  const int contID = continuums[m_continuum];
  
  xml_node<char> *node = 0;
  xml_node<char> *peak_node = doc->allocate_node( node_element, "Peak" );
  
  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  peak_node->append_attribute( att );
  
  snprintf( buffer, sizeof(buffer), "%i", contID );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "continuumID", val );
  peak_node->append_attribute( att );
  
  parent->append_node( peak_node );
  
  if( m_userLabel.size() )
  {
    val = doc->allocate_string( m_userLabel.c_str() );
    node = doc->allocate_node( node_element, "UserLabel", val );
    peak_node->append_node( node );
  }//if( m_userLabel.size() )
  
  
  switch( m_type )
  {
    case GaussianDefined: val = "GaussianDefined"; break;
    case DataDefined:     val = "DataDefined";     break;
  }//switch( m_type )
  
  node = doc->allocate_node( node_element, "Type", val );
  peak_node->append_node( node );
  
  switch( m_skewType )
  {
    case NoSkew:     val = "NoSkew";     break;
    case LandauSkew: val = "LandauSkew"; break;
  }//switch( m_skewType )
  
  node = doc->allocate_node( node_element, "Skew", val );
  peak_node->append_node( node );
  
  
  for( CoefficientType t = CoefficientType(0); 
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const char *label = to_string( t );
    
    snprintf( buffer, sizeof(buffer), "%1.8e %1.8e", m_coefficients[t], m_uncertainties[t] );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, label, val );
    
    att = doc->allocate_attribute( "fit", (m_fitFor[t] ? "true" : "false") );
    node->append_attribute( att );
    
    peak_node->append_node( node );
  }//for(...)
  
  att = doc->allocate_attribute( "forCalibration", (m_useForCalibration ? "true" : "false") );
  peak_node->append_attribute( att );
  
  att = doc->allocate_attribute( "source", (m_useForShieldingSourceFit ? "true" : "false") );
  peak_node->append_attribute( att );
  
  const char *gammaTypeVal = 0;
  switch( m_sourceGammaType )
  {
    case QLPeakDef::NormalGamma:       gammaTypeVal = "NormalGamma";       break;
    case QLPeakDef::AnnihilationGamma: gammaTypeVal = "AnnihilationGamma"; break;
    case QLPeakDef::SingleEscapeGamma: gammaTypeVal = "SingleEscapeGamma"; break;
    case QLPeakDef::DoubleEscapeGamma: gammaTypeVal = "DoubleEscapeGamma"; break;
    case QLPeakDef::XrayGamma:         gammaTypeVal = "XrayGamma";         break;
  }//switch( m_sourceGammaType )

  
  if( m_parentNuclide.size() )
  {    
    xml_node<char> *nuc_node = doc->allocate_node( node_element, "Nuclide" );
    peak_node->append_node( nuc_node );
    
    val = doc->allocate_string( m_parentNuclide.c_str() );
    node = doc->allocate_node( node_element, "Name", val );
    nuc_node->append_node( node );
    
    if( m_gammaEnergy > 0.0 )
    {
      //string transistion_parent, decay_child;

      //const SandiaDecay::Nuclide *trans_parent = m_transition->parent;
      //transistion_parent = trans_parent->symbol;
      //if( m_transition->child )
        //decay_child = m_transition->child->symbol;
      //const double energy = m_transition->products[m_radparticleIndex].energy;
      
      //val = doc->allocate_string( transistion_parent.c_str() );
      //node = doc->allocate_node( node_element, "DecayParent", val );
      //nuc_node->append_node( node );
      
      //val = doc->allocate_string( decay_child.c_str() );
      //node = doc->allocate_node( node_element, "DecayChild", val );
      //nuc_node->append_node( node );
      
      snprintf( buffer, sizeof(buffer), "%1.8e", m_gammaEnergy );
      val = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "DecayGammaEnergy", val );
      nuc_node->append_node( node );
    }//if( m_transition )
    
    node = doc->allocate_node( node_element, "DecayGammaType", gammaTypeVal );
    nuc_node->append_node( node );
  }//if( m_parentNuclide )
  
  if( m_xrayElement.size() )
  {
    xml_node<char> *xray_node = doc->allocate_node( node_element, "XRay" );
    peak_node->append_node( xray_node );
    
    val = doc->allocate_string( m_xrayElement.c_str() );
    node = doc->allocate_node( node_element, "Element", val );
    xray_node->append_node( node );
    
    snprintf( buffer, sizeof(buffer), "%1.8e", m_xrayEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "Energy", val );
    xray_node->append_node( node );
  }//if( m_xrayElement )
  
  if( m_reaction.size() )
  {
    xml_node<char> *rctn_node = doc->allocate_node( node_element, "Reaction" );
    peak_node->append_node( rctn_node );
    
    val = doc->allocate_string( m_xrayElement.c_str() );
    node = doc->allocate_node( node_element, "Name", val );
    rctn_node->append_node( node );
    
    snprintf( buffer, sizeof(buffer), "%1.8e", m_reactionEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "Energy", val );
    rctn_node->append_node( node );
    
    node = doc->allocate_node( node_element, "Type", gammaTypeVal );
    rctn_node->append_node( node );
  }//if( m_reaction )
  
  return peak_node;
}//rapidxml::xml_node<char> *toXml(...)



void QLPeakDef::fromXml( const rapidxml::xml_node<char> *peak_node,
             const std::map<int,std::shared_ptr<QLPeakContinuum> > &continuums )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !peak_node )
    throw logic_error( "QLPeakDef::fromXml(...): invalid input node" );
  
  if( !compare( peak_node->name(), peak_node->name_size(), "Peak", 4, false ) )
    throw std::logic_error( "QLPeakDef::fromXml(...): invalid input node name" );
  
  reset();
  
  int contID;
  xml_attribute<char> *att = peak_node->first_attribute( "continuumID", 11 );
  if( !att )
    throw runtime_error( "No continuum ID" );
  if( sscanf( att->value(), "%i", &contID ) != 1 )
    throw runtime_error( "Non integer continuum ID" );
  
  std::map<int,std::shared_ptr<QLPeakContinuum> >::const_iterator contpos;
  contpos = continuums.find( contID );
  if( contpos == continuums.end() )
  {
    cout << "Couldnt find valud continuum for peak: contID=" << contID << " out of " << continuums.size() << endl;
    throw runtime_error( "Couldnt find valud continuum for peak" );
  }
  
  m_continuum = contpos->second;
  
  att = peak_node->first_attribute( "forCalibration", 14 );
  if( !att )
    throw runtime_error( "missing forCalibration attribute" );
  m_useForCalibration = compare(att->value(),att->value_size(),"true",4,false);
  if( !m_useForCalibration && !compare(att->value(),att->value_size(),"false",5,false) )
    throw runtime_error( "invalis forCalibration value" );
  
  att = peak_node->first_attribute( "source", 6 );
  if( !att )
    throw runtime_error( "missing source attribute" );
  m_useForShieldingSourceFit = compare(att->value(),att->value_size(),"true",4,false);
  if( !m_useForShieldingSourceFit && !compare(att->value(),att->value_size(),"false",5,false) )
    throw runtime_error( "invalis source value" );
  
  att = peak_node->first_attribute( "version", 7 );
  if( !att )
    throw runtime_error( "missing version attribute" );
  
  int version;
  if( sscanf( att->value(), "%i", &version ) != 1 )
    throw runtime_error( "Non integer version number" );
  
  if( version != sm_xmlSerializationVersion )
    throw runtime_error( "Invalid peak version" );
  
  xml_node<char> *node = peak_node->first_node("UserLabel",9);
  if( node && node->value() )
    m_userLabel = node->value();
  
  node = peak_node->first_node("Type",4);
  if( !node || !node->value() )
    throw runtime_error( "No peak type" );
  
  if( compare(node->value(),node->value_size(),"GaussianDefined",15,false) )
    m_type = GaussianDefined;
  else if( compare(node->value(),node->value_size(),"DataDefined",11,false) )
    m_type = DataDefined;
  else
    throw runtime_error( "Invalid peak type" );
  
  node = peak_node->first_node("Skew",4);
  if( !node || !node->value() )
    throw runtime_error( "No peak skew type" );
  
  if( compare(node->value(),node->value_size(),"NoSkew",6,false) )
    m_skewType = NoSkew;
  else if( compare(node->value(),node->value_size(),"LandauSkew",10,false) )
    m_skewType = LandauSkew;
  else
    throw runtime_error( "Invalid peak skew type" );
    
  
  for( CoefficientType t = CoefficientType(0); 
      t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const char *label = to_string( t );
    
    node = peak_node->first_node(label);
    if( !node || !node->value() )
      throw runtime_error( "No coefficent " + string(label) );
    
    float dblval, dbluncrt;
    if( sscanf(node->value(), "%g %g", &dblval, &dbluncrt) != 2 )
      throw runtime_error( "unable to read value or uncert for " + string(label) );
    
    m_coefficients[t] = dblval;
    m_uncertainties[t] = dbluncrt;
    
    att = node->first_attribute("fit",3);
    if( !att || !att->value() )
      throw runtime_error( "No fit attribute for " + string(label) );
    
    m_fitFor[t] = compare(att->value(),att->value_size(),"true",4,false);
    if( !m_fitFor[t] && !compare(att->value(),att->value_size(),"false",5,false) )
      throw runtime_error( "invalid fit value" );
  }//for(...)

  xml_node<char> *nuc_node = peak_node->first_node("Nuclide",7);
  xml_node<char> *xray_node = peak_node->first_node("XRay",4);
  xml_node<char> *rctn_node = peak_node->first_node("Reaction",8);
  
  try
  {
  
    if( nuc_node )
    {
      xml_node<char> *name_node = nuc_node->first_node("Name",4);
      xml_node<char> *p_node = nuc_node->first_node("DecayParent",11);
      xml_node<char> *c_node = nuc_node->first_node("DecayChild",10);
      xml_node<char> *e_node = nuc_node->first_node("DecayGammaEnergy",16);
      xml_node<char> *type_node = nuc_node->first_node("DecayGammaType",14);
    
      //const bool isNormalNucTrans = (p_node && c_node && e_node && name_node->value()
      //                              && p_node->value() && c_node->value() && e_node->value());
    
      bool gotGammaType = false;
      const char *typeval = type_node->value();
      const size_t typelen = type_node->value_size();
    
      if( compare( typeval, typelen, "NormalGamma", 11, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = QLPeakDef::NormalGamma;
      }else if( compare( typeval, typelen, "AnnihilationGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = QLPeakDef::AnnihilationGamma;
      }else if( compare( typeval, typelen, "SingleEscapeGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = QLPeakDef::SingleEscapeGamma;
      }else if( compare( typeval, typelen, "DoubleEscapeGamma", 17, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = QLPeakDef::DoubleEscapeGamma;
      }else if( compare( typeval, typelen, "XrayGamma", 9, false ) )
      {
        gotGammaType = true;
        m_sourceGammaType = QLPeakDef::XrayGamma;
      }

      if( !name_node || !gotGammaType )
        throw runtime_error( "Invalidly specified nuclide" );
    
      m_parentNuclide = name_node->value();
      m_gammaEnergy = 0.0;
      if( e_node && sscanf( e_node->value(), "%g", &m_gammaEnergy ) != 1 )
        throw runtime_error( "Invalid nuclide gamma energy" );
    }//if( nuc_node )
    
    if( xray_node )
    {
      xml_node<char> *el_node = xray_node->first_node("Element",7);
      xml_node<char> *energy_node = xray_node->first_node("Energy",6);
  
      if( !el_node || !el_node->value() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified xray" );
    
      m_xrayElement = el_node->value();
      float dummyval;
      if( sscanf( energy_node->value(), "%g", &dummyval ) != 1 )
        throw runtime_error( "non numeric xray energy" );
      m_xrayEnergy = dummyval;
    }//if( xray_node )
  
    if( rctn_node )
    {
      xml_node<char> *name_node   = rctn_node->first_node("Name",4);
      xml_node<char> *energy_node = rctn_node->first_node("Energy",6);
      xml_node<char> *type_node   = rctn_node->first_node("Type",4);
    
      if( !name_node || !name_node->value() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified reaction" );
    
      float dummyval;
      if( sscanf( energy_node->value(), "%g", &dummyval ) != 1 )
        throw runtime_error( "non numeric reaction energy" );
      m_reactionEnergy = dummyval;
  
      //We will default do NormalGamma since early versions of serializtion didnt
      //  write the type
      m_sourceGammaType = QLPeakDef::NormalGamma;
    
      if( type_node )
      {
        const char *typeval = type_node->value();
        const size_t typelen = type_node->value_size();
  
        if( compare( typeval, typelen, "NormalGamma", 11, false ) )
          m_sourceGammaType = QLPeakDef::NormalGamma;
        else if( compare( typeval, typelen, "AnnihilationGamma", 17, false ) )
          m_sourceGammaType = QLPeakDef::AnnihilationGamma;
        else if( compare( typeval, typelen, "SingleEscapeGamma", 17, false ) )
          m_sourceGammaType = QLPeakDef::SingleEscapeGamma;
        else if( compare( typeval, typelen, "DoubleEscapeGamma", 17, false ) )
          m_sourceGammaType = QLPeakDef::DoubleEscapeGamma;
        else if( compare( typeval, typelen, "XrayGamma", 9, false ) )
          m_sourceGammaType = QLPeakDef::XrayGamma;
      }//if( type_node )

      m_reaction = name_node->value();
    }//if( rctn_node )
  }catch( std::exception &e )
  {
    m_sourceGammaType = NormalGamma;
    m_parentNuclide = m_reaction = m_xrayElement = "";
    m_gammaEnergy = m_xrayEnergy = m_reactionEnergy = 0.0;
    cerr << "Failed to assign peak at " << mean() << " keV to nuclide/xray/reaction: " << e.what();
  }
}//void fromXml(...)



std::shared_ptr<QLPeakContinuum> QLPeakDef::continuum()
{
  return m_continuum;
}

std::shared_ptr<const QLPeakContinuum> QLPeakDef::continuum() const
{
  return m_continuum;
}

std::shared_ptr<QLPeakContinuum> QLPeakDef::getContinuum()
{
  return m_continuum;
}

void QLPeakDef::setContinuum( std::shared_ptr<QLPeakContinuum> continuum )
{
  if( !continuum )
    throw runtime_error( "QLPeakDef::setContinuum(...): invalid input" );
  m_continuum = continuum;
}//void setContinuum(...)


void QLPeakDef::makeUniqueNewContinuum()
{
  m_continuum = std::make_shared<QLPeakContinuum>(*m_continuum);
}


//The bellow should in principle take care of gaussian area and the skew area
double QLPeakDef::peakArea() const
{
  double amp = m_coefficients[QLPeakDef::GaussAmplitude];
  
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
    break;
    
    case QLPeakDef::LandauSkew:
      amp += skew_integral( lowerX(), upperX() );
    break;
  }//switch( m_skewType )
  
  return amp;
}//double peakArea() const


double QLPeakDef::peakAreaUncert() const
{
  double uncert = m_uncertainties[QLPeakDef::GaussAmplitude];
  
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
    break;
      
    case QLPeakDef::LandauSkew:
    {
      const double skew_area = skew_integral( lowerX(), upperX() );
      const double frac_uncert = m_uncertainties[QLPeakDef::LandauAmplitude]
                                 / m_coefficients[QLPeakDef::LandauAmplitude];
      const double skew_uncert = skew_area * frac_uncert;
      //XXX - bellow assumes ampltude and skew amplitude are uncorelated,
      //      which they are not.
      uncert = sqrt( uncert*uncert + skew_uncert*skew_uncert );
      break;
    }//case QLPeakDef::LandauSkew:
  }//switch( m_skewType )
  
  return uncert;
}//double peakAreaUncert() const


void QLPeakDef::setPeakArea( const double a )
{
  double area = a; //m_coefficients[QLPeakDef::GaussAmplitude];
  
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
    break;
      
    case QLPeakDef::LandauSkew:
    {
      const double skew_area = skew_integral( lowerX(), upperX() );
      const double skew_frac = skew_area / (skew_area + area);
      area = (1.0 - skew_frac) * a;
      break;
    }//case QLPeakDef::LandauSkew:
  }//switch( m_skewType )
  
  m_coefficients[QLPeakDef::GaussAmplitude] = area;
}//void QLPeakDef::setPeakArea( const double a )


void QLPeakDef::setPeakAreaUncert( const double uncert )
{
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
      m_uncertainties[QLPeakDef::GaussAmplitude] = uncert;
    break;
      
    case QLPeakDef::LandauSkew:
    {
      //XXX - the bellow only modifies the gaus uncertainty, and not the skew
      //  uncertainty - I was just to lazy to do it properly (should also look
      //  at how coorelated the skew amplitude error is to peak amplitude error
      //  ...)
      const double skew_area = skew_integral( lowerX(), upperX() );
      const double gauss_area = m_coefficients[QLPeakDef::GaussAmplitude];
      const double total_area = skew_area + gauss_area;
      const double frac_uncert = uncert / total_area;
      const double new_gauss_uncert = frac_uncert * gauss_area;
      m_uncertainties[QLPeakDef::GaussAmplitude] = new_gauss_uncert;
      
      /*
      const double gaussuncert = m_coefficients[QLPeakDef::GaussAmplitude];
      const double skewuncert = m_coefficients[QLPeakDef::LandauAmplitude];
      
      if( skewuncert > 0.0 && skew_area > 0.0 )
      {
        const double gauss_area = m_coefficients[QLPeakDef::GaussAmplitude];
        const double area = skew_area + gauss_area;
        const double skew_uncert = skew_area * skewuncert
                                   / m_coefficients[QLPeakDef::LandauAmplitude];
        const double skew_frac_uncert = skew_uncert / (skew_uncert+gaussuncert);
        
        const double skew_uncert = skew_area * skew_frac_uncert
        
        m_uncertainties[QLPeakDef::GaussAmplitude] = fracuncert*gauss_area;
        m_uncertainties[QLPeakDef::LandauAmplitude] = fracuncert*;
      }else
      {
        m_uncertainties[QLPeakDef::GaussAmplitude] = uncert;
      }//if( skewuncert > 0.0 ) / else
       */
      break;
    }//case QLPeakDef::LandauSkew:
  }//switch( m_skewType )
}//void setPeakAreaUncert( const double a )



bool QLPeakDef::lessThanByMean( const QLPeakDef &lhs, const QLPeakDef &rhs )
{
  return (lhs.m_coefficients[QLPeakDef::Mean] < rhs.m_coefficients[QLPeakDef::Mean]);
}//lessThanByMean(...)

bool QLPeakDef::lessThanByMeanShrdPtr( const std::shared_ptr<const QLPeakDef> &lhs,
                                  const std::shared_ptr<const QLPeakDef> &rhs )
{
  if( !lhs || !rhs )
    return (lhs < rhs);
  return lessThanByMean( *lhs, *rhs );
}




bool QLPeakDef::operator==( const QLPeakDef &rhs ) const
{
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    if( m_coefficients[t] != rhs.m_coefficients[t] )
      return false;
  }

  return m_type==rhs.m_type
         && (*m_continuum == *rhs.m_continuum)
      && m_parentNuclide==rhs.m_parentNuclide
      && m_gammaEnergy == rhs.m_gammaEnergy
      && m_sourceGammaType==rhs.m_sourceGammaType
      && m_xrayElement==rhs.m_xrayElement
      && m_xrayEnergy==rhs.m_xrayEnergy
      && m_reaction==rhs.m_reaction
      && m_reactionEnergy==rhs.m_reactionEnergy
      && m_useForCalibration==rhs.m_useForCalibration
      && m_useForShieldingSourceFit==rhs.m_useForShieldingSourceFit
      ;
}//QLPeakDef::operator==


std::string QLPeakDef::parentNuclide() const
{
  return m_parentNuclide;
}//string parentNuclide() const;


QLPeakDef::SourceGammaType QLPeakDef::sourceGammaType() const
{
  return m_sourceGammaType;
}

std::string QLPeakDef::xrayElement() const
{
  return m_xrayElement;
}

float QLPeakDef::xrayEnergy() const
{
  return m_xrayEnergy;
}

std::string QLPeakDef::reaction() const
{
  return m_reaction;
}

float QLPeakDef::reactionEnergy() const
{
  return m_reactionEnergy;
}


float QLPeakDef::gammaParticleEnergy() const
{
  if( m_parentNuclide.size() )
  {
    switch( m_sourceGammaType )
    {
      case NormalGamma:
      case XrayGamma:
        return m_gammaEnergy;
      break;
      case AnnihilationGamma:
        return 510.99891f;
      case SingleEscapeGamma:
        return m_gammaEnergy - 510.9989f;
      break;
      case DoubleEscapeGamma:
        return m_gammaEnergy - 2.0f*510.9989f;
      break;
    }
  }//if( m_parentNuclide )
  
  if( m_xrayElement.size() )
    return m_xrayEnergy;
  
  if( m_reaction.size() )
    return m_reactionEnergy;
  
  //throw runtime_error( "Peak doesnt have a gamma associated with it" );
  
  return 0.0f;
}//float gammaParticleEnergy() const


bool QLPeakContinuum::defined() const
{
  switch( m_type )
  {
    case NoOffset:
      return false;
      
    case Cubic:
      if( m_values[3] != 0.0)
        return true;
    case Quardratic:
      if( m_values[2] != 0.0)
        return true;
    case Linear:
      if( m_values[1] != 0.0)
        return true;
    case Constant:
      return (m_values[0] != 0.0);
    break;
    
    case External:
      return !!m_externalContinuum;
    break;
  }//switch( m_type )
  
  return false;
}//bool defined() const




double QLPeakDef::lowerX() const
{
  if( m_continuum->QLPeakContinuum::energyRangeDefined() )
    return m_continuum->lowerEnergy();

  if( m_skewType == QLPeakDef::LandauSkew )
    return min(m_coefficients[QLPeakDef::Mean] - m_coefficients[QLPeakDef::LandauSigma]+(0.22278*m_coefficients[QLPeakDef::LandauSigma])-25.0*m_coefficients[QLPeakDef::LandauSigma],
        m_coefficients[QLPeakDef::Mean] - 4.0*m_coefficients[QLPeakDef::Sigma]);

  return m_coefficients[QLPeakDef::Mean] - 4.0*m_coefficients[QLPeakDef::Sigma];
}//double lowerX() const


double QLPeakDef::upperX() const
{
  if( m_continuum->QLPeakContinuum::energyRangeDefined() )
    return m_continuum->upperEnergy();
  
  return m_coefficients[QLPeakDef::Mean] + 4.0*m_coefficients[QLPeakDef::Sigma];
}//double upperX() const


double QLPeakDef::gauss_integral( const double x0, const double x1 ) const
{
  double integral = gaus_integral( m_coefficients[QLPeakDef::Mean],
                                   m_coefficients[QLPeakDef::Sigma],
                                   m_coefficients[QLPeakDef::GaussAmplitude],
                                   x0, x1 );
  
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
    break;
    
    case QLPeakDef::LandauSkew:
      integral += skew_integral( x0, x1 );
    break;
  };//enum SkewType

  return integral;
}//double gauss_integral( const double x0, const double x1 ) const;


double QLPeakDef::offset_integral( const double x0, const double x1 ) const
{
  return m_continuum->offset_integral( x0, x1 );
}//double offset_integral( const double x0, const double x1 ) const


double QLPeakDef::landau_potential_lowerX( const double peak_mean,
                                          const double peak_sigma )
{
  return peak_mean - 8.0*peak_sigma;
}//double landau_potential_lowerX() const


double QLPeakDef::landau_potential_upperX( const double peak_mean,
                                          const double /*peak_sigma*/)
{
  return peak_mean;
}//double landau_potential_upperX() const


double QLPeakDef::landau_integral( const double x0, const double x1,
                                  const double peak_mean,
                                  const double amplitude,
                                  const double mode,
                                  const double sigma )
{
  if( amplitude <= 0.0 )
    return 0.0;

  const double y0 = landau_cdf( peak_mean - x0, mode, sigma );
  const double y1 = landau_cdf( peak_mean - x1, mode, sigma );
  
  return amplitude*(y0-y1);
}//static double landau_integral( ... )

double QLPeakDef::skew_integral( const double xbinlow, const double xbinup,
                               const double peak_amplitude,
                               const double peak_mean,
                               const double t, const double s, const double b )
{
  //XXX - should make landaumode and landausigma relative to peak sigma
  return peak_amplitude*landau_integral( xbinlow, xbinup, peak_mean, t, s, b );

  /*
  //Below is an implementation of the skew definition TSpectrumFit uses.
  //  I cant seem to get it to work real well when using TMinuit2 to fit for
  //  the skew - additionally I dont know if the s*erfc( p )/2 term is well
  //  motoivated.
  const double amp = peak_amplitude /sqrt( 2.0*boost::math::constants::pi<double>()*sigma*sigma );

  const double x = (xbinup+xbinlow)/2.0;
  const double p = (x-peak_mean)/sigma;
  const double c = p + 1.0 / (2.0 * b);
  double e = p / b;
  if( e > 700.0 )
    e = 700.0;

  double skew = 0.0;
  if( b != 0.0 )
    skew += 0.5*t*exp( e ) * boost::math::erfc( c );
  skew += 0.5*s*boost::math::erfc( p );

  return amp * skew;
  */
}//double QLPeakDef::skew_integral(...)


double QLPeakDef::skew_integral( const double x0, const double x1 ) const
{
  double area = 0.0;
  
  switch( m_skewType )
  {
    case QLPeakDef::NoSkew:
      break;
      
    case QLPeakDef::LandauSkew:
      area = skew_integral( x0, x1, m_coefficients[QLPeakDef::GaussAmplitude],
                           m_coefficients[QLPeakDef::Mean],
                           m_coefficients[QLPeakDef::LandauAmplitude],
                           m_coefficients[QLPeakDef::LandauMode],
                           m_coefficients[QLPeakDef::LandauSigma] );
      break;
  };//switch( m_skewType )
  
  return area;
}//double skew_integral( const double x0, const double x1 ) const



double QLPeakDef::gaus_integral( const double peak_mean, const double peak_sigma,
                               const double peak_amplitude,
                               const double x0, const double x1 )
{
  if( peak_sigma==0.0 || peak_amplitude==0.0 )
    return 0.0;

  //Since the data is only float accuracy, there is no reason to calculate
  //  things past 9 digits (the erf function does hit the profiler pretty
  //  decently while fitting for peaks)
  using boost::math::policies::digits10;
  typedef boost::math::policies::policy<digits10<9> > my_pol_9;
  
  const double sqrt2 = boost::math::constants::root_two<double>();
  
  const double range = x1 - x0;
  const double erflowarg = (x0-peak_mean)/(sqrt2*peak_sigma);
  const double erfhigharg = (x0+range-peak_mean)/(sqrt2*peak_sigma);
  const double cdflow = 0.5*( 1.0 + boost::math::erf( erflowarg, my_pol_9() ) );
  const double cdhigh = 0.5*( 1.0 + boost::math::erf( erfhigharg, my_pol_9() ) );

  return peak_amplitude * (cdhigh - cdflow);
}//double gaus_integral(...)


////////////////////////////////////////////////////////////////////////////////

QLPeakContinuum::QLPeakContinuum()
: m_type( QLPeakContinuum::NoOffset ),
  m_lowerEnergy( 0.0 ),
  m_upperEnergy( 0.0 ),
  m_refernceEnergy( 0.0 )
{
}//QLPeakContinuum constructor

bool QLPeakContinuum::operator==( const QLPeakContinuum &rhs ) const
{
  return m_type==rhs.m_type
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_upperEnergy == rhs.m_upperEnergy
         && m_refernceEnergy == rhs.m_refernceEnergy
         && m_values == rhs.m_values
         && m_uncertainties == rhs.m_uncertainties
         && m_fitForValue == rhs.m_fitForValue
         && m_externalContinuum == rhs.m_externalContinuum;
}

void QLPeakContinuum::setParameters( double referenceEnergy,
                                   const std::vector<double> &x,
                                   const std::vector<double> &uncertainties )
{
  switch( m_type )
  {
    case NoOffset: case External:
      throw runtime_error( "QLPeakContinuum::setParameters invalid m_type" );
      
    case Constant:   case Linear:
    case Quardratic: case Cubic:
    break;
  };//switch( m_type )
  
  if( x.size() != static_cast<size_t>(m_type) )
    throw runtime_error( "QLPeakContinuum::setParameters invalid parameter size" );
  
  m_values = x;
  m_refernceEnergy = referenceEnergy;
  m_fitForValue.resize( m_values.size(), true );
  
  if( uncertainties.empty() )
  {
    m_uncertainties.clear();
    m_uncertainties.resize( x.size(), 0.0 );
  }else
  {
    if( uncertainties.size() != static_cast<size_t>(m_type) )
      throw runtime_error( "QLPeakContinuum::setParameters invalid uncert size" );
    m_uncertainties = uncertainties;
  }//if( uncertainties.empty() ) / else
}//void setParameters(...)



bool QLPeakContinuum::setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )
{
  if( polyCoefNum >= m_fitForValue.size() )
    return false;
  
  m_fitForValue[polyCoefNum] = fit;
  return true;
}//bool setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )


bool QLPeakContinuum::setPolynomialCoef( size_t polyCoef, double val )
{
  if( polyCoef >= m_values.size() )
    return false;
  
  m_values[polyCoef] = val;
  return true;
}


bool QLPeakContinuum::setPolynomialUncert( size_t polyCoef, double val )
{
  if( polyCoef >= m_uncertainties.size() )
    return false;
  
  m_uncertainties[polyCoef] = val;
  return true;
}

void QLPeakContinuum::setParameters( double referenceEnergy,
                                   const double *parameters,
                                   const double *uncertainties )
{
  switch( m_type )
  {
    case NoOffset: case External:
      throw runtime_error( "QLPeakContinuum::setParameters invalid m_type" );
      
    case Constant:   case Linear:
    case Quardratic: case Cubic:
      break;
  };//switch( m_type )
  
  if( !parameters )
    throw runtime_error( "QLPeakContinuum::setParameters invalid paramters" );
  
  vector<double> uncerts, values( parameters, parameters+m_type );
  if( uncertainties )
    uncerts.insert( uncerts.end(), uncertainties, uncertainties+m_type );
  
  setParameters( referenceEnergy, values, uncerts );
}//setParameters


void QLPeakContinuum::setExternalContinuum( const std::shared_ptr<const Measurement> &data )
{
  if( m_type != External )
    throw runtime_error( "QLPeakContinuum::setExternalContinuum invalid m_type" );
 
  m_externalContinuum = data;
}//setExternalContinuum(...)


void QLPeakContinuum::setRange( const double lower, const double upper )
{
  m_lowerEnergy = lower;
  m_upperEnergy = upper;
  if( m_lowerEnergy > m_upperEnergy )
    std::swap( m_lowerEnergy, m_upperEnergy );
}//void setRange( const double lowerenergy, const double upperenergy )


bool QLPeakContinuum::energyRangeDefined() const
{
  return (m_lowerEnergy != m_upperEnergy);
}//bool energyRangeDefined() const


bool QLPeakContinuum::isPolynomial() const
{
  switch( m_type )
  {
    case NoOffset:   case External:
      return false;
    case Constant:   case Linear:
    case Quardratic: case Cubic:
      return true;
    break;
  }//switch( m_type )
  
  return false;
}//bool isPolynomial() const


void QLPeakContinuum::setType( QLPeakContinuum::OffsetType type )
{
  m_type = type;
  
  switch( m_type )
  {
    case NoOffset:
      m_values.clear();
      m_uncertainties.clear();
      m_fitForValue.clear();
      m_externalContinuum.reset();
      m_lowerEnergy = m_upperEnergy = m_refernceEnergy = 0.0;
    break;
      
    case Constant:
      m_values.resize( 1, 0.0 );
      m_uncertainties.resize( 1, 0.0 );
      m_fitForValue.resize( 1, true );
      m_externalContinuum.reset();
    break;
    
    case Linear:
      m_values.resize( 2, 0.0 );
      m_uncertainties.resize( 2, 0.0 );
      m_fitForValue.resize( 2, true );
      m_externalContinuum.reset();
    break;
      
    case Quardratic:
      m_values.resize( 3, 0.0 );
      m_uncertainties.resize( 3, 0.0 );
      m_fitForValue.resize( 3, true );
      m_externalContinuum.reset();
    break;
    
    case Cubic:
      m_values.resize( 4, 0.0 );
      m_uncertainties.resize( 4, 0.0 );
      m_fitForValue.resize( 4, true );
      m_externalContinuum.reset();
    break;
      
    case External:
      m_values.clear();
      m_uncertainties.clear();
      m_fitForValue.clear();
      m_refernceEnergy = 0.0;
    break;
  };//switch( m_type )
  
}//void setType( QLPeakContinuum::OffsetType type )

void QLPeakContinuum::calc_linear_continuum_eqn( const std::shared_ptr<const Measurement> &data,
                                        const double x0, const double x1,
                                        const int nbin )
{
  m_refernceEnergy = x0;
  m_lowerEnergy = x0;
  m_upperEnergy = x1;
  m_type = QLPeakContinuum::Linear;
  m_values.resize( 2 );
  m_fitForValue.resize( 2, true );
  
  const size_t xlowchannel = data->find_gamma_channel( (float)x0 );
  const size_t xhighchannel = data->find_gamma_channel( (float)x1 );
  
  eqn_from_offsets( xlowchannel, xhighchannel, m_refernceEnergy, data, nbin, m_values[1], m_values[0] );
}//void calc_continuum_eqn(...)


double QLPeakContinuum::offset_integral( const double x0, const double x1 ) const
{
  switch( m_type )
  {
    case NoOffset:
      return 0.0;
    case Constant:   case Linear:
    case Quardratic: case Cubic:
      return offset_eqn_integral( &(m_values[0]),
                                 m_type, x0, x1, m_refernceEnergy );
    case External:
      if( !m_externalContinuum )
        return 0.0;
      return gamma_integral( m_externalContinuum, x0, x1 );
    break;
  };//enum OffsetType
  
  return 0.0;
}//double offset_integral( const double x0, const double x1 ) const


void QLPeakContinuum::eqn_from_offsets( size_t lowchannel,
                                      size_t highchannel,
                               const double peakMean,
                               const std::shared_ptr<const Measurement> &data,
                               const size_t nbin,
                               double &m, double &b )
{
  double x1  = data->gamma_channel_lower( lowchannel ) - peakMean;
  const double dx1 = data->gamma_channel_width( lowchannel );
  
  //XXX - hack! Bellow 'c' will be inf if the following check would fail; I
  //      really should work out how to fix this properly with math, but for
  //      right now I'll fudge it, which will toss the answer off a bit, but
  //      whatever.
  if( fabs(2.0*x1*dx1 + dx1*dx1) < DBL_EPSILON )
  {
    static int ntimes = 0;
    if( ntimes++ < 4 )
      cerr << "Should fix math in QLPeakContinuum::eqn_from_offsets" << endl;
    x1 = data->gamma_channel_lower( lowchannel+1 ) - peakMean;
  }
  
  const double x2  = data->gamma_channel_lower( highchannel ) - peakMean;
  
  const double dx2 = data->gamma_channel_width( highchannel );
  
  if( lowchannel < nbin )
    lowchannel = nbin;
  highchannel = std::min( highchannel, data->num_gamma_channels()-1 );
  
  const double nbinInv = 1.0 / (1.0+2.0*nbin);
  const double y1 = nbinInv*data->gamma_channels_sum( lowchannel-nbin, lowchannel+nbin );
  const double y2 = nbinInv*data->gamma_channels_sum( highchannel-nbin, highchannel+nbin);

  const double c = (2.0*x2*dx2 + dx2*dx2)/(2.0*x1*dx1 + dx1*dx1);
  b = (y2-y1*c)/(dx2-dx1*c);
  m = 2.0*(y1-b*dx1)/(2.0*x1*dx1+dx1*dx1);
  
  if( IsNan(m) || IsInf(m) || IsNan(b) || IsInf(b) )
  {
    cerr << "QLPeakContinuum::eqn_from_offsets(...): Invalid results" << endl;
    m = b = 0.0;
  }
}//eqn_from_offsets(...)


double QLPeakContinuum::offset_eqn_integral( const double *coefs,
                                           QLPeakContinuum::OffsetType type,
                                           double x0, double x1,
                                           const double peak_mean )
{
  switch( type )
  {
    case NoOffset: case External:
      throw runtime_error( "QLPeakContinuum::offset_eqn_integral(...) may only be"
                           " called for polynomial backgrounds" );
      
    case Constant:   case Linear:
    case Quardratic: case Cubic:
    break;
  };//enum OffsetType
    
  x0 -= peak_mean;
  x1 -= peak_mean;
  
  double answer = 0.0;
  const int maxorder = static_cast<int>( type );
  for( int order = 0; order < maxorder; ++order )
  {
    const double exp = order + 1.0;
    answer += (coefs[order]/exp) * (std::pow(x1,exp) - std::pow(x0,exp));
  }//for( int order = 0; order < maxorder; ++order )
  
  return std::max( answer, 0.0 );
}//offset_eqn_integral(...)


void QLPeakContinuum::translate_offset_polynomial( double *new_coefs,
                                                 const double *old_coefs,
                                                 QLPeakContinuum::OffsetType type,
                                                 const double new_center,
                                                 const double old_center )
{
  switch( type )
  {
    case NoOffset:
    case External:
      throw runtime_error( "translate_offset_polynomial invalid offset type" );
      
    case Quardratic:
    case Cubic:
      throw runtime_error( "translate_offset_polynomial does not yet support "
                           "quadratic or cubic polynomials" );
      
    case Linear:
    {
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      break;
    }//case Linear:
      
    case Constant:
      new_coefs[0] = old_coefs[0];
      break;
  }//switch( type )
}//void translate_offset_polynomial(...)

